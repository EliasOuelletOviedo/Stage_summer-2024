"""Étape 1 — décodage : chaque photon dans sa passe, son genre et son canal de temps.

Usage :  python3 step1_decode.py <dossier>          (~1 min)

Une passe = un scan de 210 ms (créneau de 240 ms de la NI), bornée par les M0/M3 de la
carte. Dans chaque passe, le code de routage sépare les échantillons d'acquisition (code
1 : laser opto éteint, photons analysés par FLIMApp) des phases opto du PWM (code 0 : le
laser à 1064 nm à sa commande, photons jetés par FLIMApp) ; les photons des pauses (30 ms,
code 0) sont gardés à part. Les stims viennent de la commande 1 relue par la NI
(readback.bin, 10 kHz) : une passe est « stim » si le laser opto est commandé pendant son
scan. Le créneau s de la NI est la passe s + 1 de la carte.

Contrôles :
  - les photons de code 1 de chaque passe doivent être exactement ceux de frames.csv ;
  - l'attribution des canaux (config.CANAL) : le CFP doit avoir le temps moyen le plus court ;
  - les galvos (readback) : immobiles si le mode ROI était éteint.

La ROI de chaque passe est le code de routage de ses photons d'acquisition (en mode ROI,
l'indice de la ROI visitée ; 1 = `CODE_SANS_ROI` quand le mode ROI est éteint) : l'étape 6
sépare les ROIs avec.

Sorties : mesures/passes.npz, results/stims.csv, qc/decodage.json, qc/validation_flimapp.csv
"""
import json
import time

import numpy as np
import pandas as pd

import config as C
import spc as S

N_PHASE = 200           # tranches de 1 µs dans la période du PWM (2 échantillons de 100 µs)


def commandes_par_creneau(r, noms, ses, n):
    """Par créneau de la NI : part des échantillons du scan où la commande 1 dépasse 0,5 V,
    position moyenne et écart-type des galvos pendant le scan."""
    c1, gx, gy = (r[:, noms.index(k)] for k in ("command_1", "galvo_x", "galvo_y"))
    lent = ses["scan"] + ses["pause"]
    entree = ses["pause"]                       # l'entrée sur la première ROI : une pause
    out = []
    for s in range(n):
        a = entree + s * lent
        b = a + ses["scan"]
        out.append(dict(creneau=s, opto_frac=float((c1[a:b] > 0.5).mean()), commande_max_v=float(c1[a:b].max()),
                        opto_pause_frac=float((c1[b:a + lent] > 0.5).mean()),
                        galvo_x_v=float(gx[a:b].mean()), galvo_y_v=float(gy[a:b].mean()),
                        galvo_sd_v=float(np.hypot(gx[a:b].std(), gy[a:b].std()))))
    return pd.DataFrame(out)


def main():
    root = C.racine()
    D = C.dirs(root)
    ses = C.session(root)
    t0 = time.time()
    r, noms = C.readback(root)
    fr = C.frames(root)
    vi = C.visits(root)
    n_creneaux = len(vi)
    K = commandes_par_creneau(r, noms, ses, n_creneaux)
    print(f"  readback : {len(r)} échantillons ({len(r) / ses['fs']:.1f} s), {n_creneaux} créneaux ; "
          f"galvos {K.galvo_x_v.mean() * 1e3:+.2f} / {K.galvo_y_v.mean() * 1e3:+.2f} mV, "
          f"écart-type {K.galvo_sd_v.max() * 1e3:.2f} mV au plus ; mode ROI {'actif' if ses['rois_actives'] else 'éteint'}")
    stim_creneau = (K.opto_frac > 0).values
    if not np.array_equal(stim_creneau, (vi.command1_v > 0).values):
        print("  ATTENTION : commande 1 relue ≠ visits.csv")

    res, bornes_ref = {}, None
    for nom, k in C.CANAL.items():
        p = ses["spc"][k]["spc"]
        m0, m3, n_ph, n_gap = S.marqueurs(p)
        debut, fin = S.bornes(m0, m3)
        if bornes_ref is None:
            bornes_ref = (debut, fin)
        elif not (np.array_equal(debut, bornes_ref[0]) and np.array_equal(fin, bornes_ref[1])):
            raise SystemExit("les deux cartes n'ont pas les mêmes passes")
        if len(debut) != n_creneaux:
            raise SystemExit(f"{len(debut)} passes pour {n_creneaux} créneaux")
        tics_ech = float(np.median(np.diff(debut))) / (ses["scan"] + ses["pause"])
        H, prof, codes, autres = S.ranger(p, debut, fin, stim_creneau, tics_ech, n_phase=N_PHASE)
        res[nom] = dict(H=H, prof=prof, codes=codes, autres=autres, n_ph=n_ph, n_gap=n_gap, tics_ech=tics_ech,
                        m0=len(m0), m3=len(m3))
        print(f"  {nom.upper()} (IN{k}) : {n_ph} photons, {len(debut)} passes, GAP {n_gap} ; "
              f"dans les passes : acquisition {H[:, S.ACQ].sum()}, phases opto {H[:, S.OPTO].sum()} ; "
              f"pauses {H[:, S.PAUSE].sum()} ; {autres}  ({time.time() - t0:.0f} s)")

    debut, fin = bornes_ref
    tic = ses["tic_s"]
    # --- contrôle : frames.csv
    val = pd.DataFrame({"passe": fr["pass"]})
    ok_all = True
    for nom, k in C.CANAL.items():
        n_moi = res[nom]["H"][:, S.ACQ].sum(1)
        val[f"photons_{nom}"] = n_moi
        val[f"flimapp_{nom}"] = fr[f"photons_ch{k}"].values
        ok = np.array_equal(n_moi, fr[f"photons_ch{k}"].values.astype(int))
        ok_all &= ok
        print(f"  contrôle frames.csv, {nom.upper()} : {'identique' if ok else 'DIFFÉRENT'} "
              f"({int((n_moi != fr[f'photons_ch{k}'].values).sum())} passes diffèrent)")
    dt_debut = np.abs(debut * tic - fr.pass_start_s.values).max()
    val.to_csv(D["qc"] / "validation_flimapp.csv", index=False)
    # --- contrôle : attribution des canaux
    t = C.temps_canaux()
    tm = {nom: float((t * res[nom]["H"][:, S.ACQ].sum(0)).sum() / res[nom]["H"][:, S.ACQ].sum()) for nom in C.CANAL}
    if tm["cfp"] >= tm["yfp"]:
        raise SystemExit(f"temps moyens CFP {tm['cfp']:.3f} ≥ YFP {tm['yfp']:.3f} ns : vérifier config.CANAL")
    print(f"  temps moyen d'arrivée : CFP {tm['cfp']:.3f} ns, YFP {tm['yfp']:.3f} ns ; "
          f"rapport des photons YFP/CFP {res['yfp']['H'][:, S.ACQ].sum() / res['cfp']['H'][:, S.ACQ].sum():.3f}")

    # --- ROI de chaque passe : le code de routage majoritaire de ses photons d'acquisition (0 : aucun)
    codes = res["cfp"]["codes"] + res["yfp"]["codes"]
    roi_passe = np.where(codes[:, 1:].sum(1) > 0, np.argmax(codes[:, 1:], 1) + 1, 0)
    rois = {int(r): int((roi_passe == r).sum()) for r in np.unique(roi_passe)}
    melange = float(1 - codes[:, 1:].max(1).sum() / max(codes[:, 1:].sum(), 1))
    print(f"  ROIs (codes de routage) : {rois} ; photons d'un autre code que celui de leur passe : {100 * melange:.3f} %")

    # --- stims : blocs de passes consécutives
    dd = np.diff(np.r_[0, stim_creneau.astype(int), 0])
    a_s, b_s = np.flatnonzero(dd == 1), np.flatnonzero(dd == -1) - 1
    stims = pd.DataFrame(dict(stim=np.arange(1, len(a_s) + 1), passe_debut=a_s + 1, passe_fin=b_s + 1,
                              n_passes=b_s - a_s + 1, t_debut_s=debut[a_s] * tic, t_fin_s=fin[b_s] * tic,
                              commande_v=[float(K.commande_max_v[a:b + 1].max()) for a, b in zip(a_s, b_s)],
                              part_opto=[float(K.opto_frac[a:b + 1].mean()) for a, b in zip(a_s, b_s)]))
    stims.to_csv(D["results"] / "stims.csv", index=False)
    print(stims.round(3).to_string(index=False))

    np.savez_compressed(
        D["mesures"] / "passes.npz",
        H_cfp=res["cfp"]["H"].astype(np.int32), H_yfp=res["yfp"]["H"].astype(np.int32),
        prof_cfp=res["cfp"]["prof"], prof_yfp=res["yfp"]["prof"],
        debut=debut, fin=fin, t_debut_s=debut * tic, t_fin_s=fin * tic, t_s=(debut + fin) / 2 * tic,
        stim=stim_creneau, opto_frac=K.opto_frac.values, tic_s=tic, roi_passe=roi_passe,
        tics_echantillon=res["cfp"]["tics_ech"], n_phase=N_PHASE,
        pwm=np.array(ses["pwm"]), scan=ses["scan"], pause=ses["pause"], fs=ses["fs"],
        flimapp_tau_cfp=fr[f"lifetime_ch{C.CANAL['cfp']}_ns"].values,
        flimapp_tau_yfp=fr[f"lifetime_ch{C.CANAL['yfp']}_ns"].values)
    K.to_csv(D["qc"] / "commandes_creneaux.csv", index=False)
    info = dict(
        passes=len(debut), validation_frames_csv=bool(ok_all), ecart_max_debut_s=float(dt_debut),
        tics_par_echantillon=res["cfp"]["tics_ech"], echantillon_us=res["cfp"]["tics_ech"] * tic * 1e6,
        mode_roi_actif=ses["rois_actives"], rois_passes=rois, part_photons_autre_code=melange,
        galvos_mV=[float(K.galvo_x_v.mean() * 1e3), float(K.galvo_y_v.mean() * 1e3)],
        galvo_sd_max_mV=float(K.galvo_sd_v.max() * 1e3), temps_moyen_ns=tm,
        **{f"{nom}_{k}": (int(v) if isinstance(v, (int, np.integer)) else v)
           for nom in C.CANAL for k, v in [("photons", res[nom]["n_ph"]), ("gap", res[nom]["n_gap"]),
                                           ("m0", res[nom]["m0"]), ("m3", res[nom]["m3"])]},
        **{f"{nom}_{k}": int(v) for nom in C.CANAL for k, v in res[nom]["autres"].items()})
    json.dump(info, open(D["qc"] / "decodage.json", "w"), indent=1, ensure_ascii=False)
    print(f"  -> mesures/passes.npz, results/stims.csv, qc/decodage.json ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
