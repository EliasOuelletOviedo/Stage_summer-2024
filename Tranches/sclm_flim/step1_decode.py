"""Étape 1 — décodage des .spc, fenêtres de stim, dérive, recalage des images Bliq.

Usage :  python3 step1_decode.py <dossier du jeu de données>

1. Chaque .spc est relu photon par photon (spc.Rangement, le rangement de FLIMApp) :
   image d'intensité et temps moyen par pixel, photons par image (16,7 Hz), et par
   bloc de `images_par_bloc` images (le bloc de FLIMApp, 50 images = 3 s) une image
   réduite 4 x 4 et le déclin de tout le champ.
2. Les fenêtres de commande du protocole (visits.csv, horloge NI) sont recalées sur
   l'horloge des cartes par la fuite de la stim dans ch2, image par image.
3. La dérive du champ est mesurée bloc par bloc (corrélation de phase).
4. Chaque image fixe Bliq est recalée sur l'image FLIM : échelle en x (pixel de 51,8 ns
   contre 51,2 ns) puis translation, au meilleur NCC.

Sorties : projections/flim.npz, projections/fovs.npz, projections/session.json,
projections/images.png, qc/par_image.png, qc/derive.png, qc/recalage_fovs.png.
"""
import json
import time

import numpy as np
import tifffile
from scipy import ndimage as ndi

import config
import figures
import spc
from recalage import _prep, correlation_phase, ncc
from config import BIN_DERIVE, CANAUX, ECHELLES_X, MIN_NCC_MEME_PLAN
from figures import INK, INK_2, SURFACE, bandes_stim, contraste, plt, sauver, style
from flim_fit import Modele

MAX_TRAMES = 400_000          # 6,6 h à 16,7 Hz
# Recalage NI -> cartes par la fuite de la stim : décalage cherché à ±RECALAGE_MAX_S (les
# sessions où la fuite est nette donnent +50 à +140 ms). Sous MIN_CORR_FUITE, la fuite est
# trop faible pour dire quoi que ce soit (stim à 25 % : corrélation 0,05) : l'horloge de la
# NI est gardée telle quelle, à ~0,1 s près.
RECALAGE_MAX_S = 0.5
MIN_CORR_FUITE = 0.1


# ------------------------------------------------------------------ décodage
def decoder_canal(info, fpb):
    ini = spc.lire_ini(info["ini"])
    geo = spc.geometrie(info["ini"])
    nb = int(ini["clamp"]["blocs"])                 # blocs complets publiés par FLIMApp
    ny, nx = geo.lignes_image, geo.pixels_par_ligne
    by, bx = ny // BIN_DERIVE, nx // BIN_DERIVE
    inten = np.zeros(ny * nx, np.int64)
    somme_c = np.zeros(ny * nx, np.float64)
    par_trame = np.zeros(MAX_TRAMES, np.int64)
    blocs = np.zeros(nb * by * bx, np.int64)
    declin = np.zeros(nb * config.CANAUX_TEMPS, np.int64)
    trames_t = []
    t0 = time.time()
    for ph in spc.Rangement(info["spc"], geo, canaux=config.CANAUX_TEMPS):
        trames_t.append(ph["trames_t"])
        tr, y, x, c = ph["trame"], ph["y"], ph["x"], ph["canal"]
        pix = y * nx + x
        inten += np.bincount(pix, minlength=ny * nx)
        somme_c += np.bincount(pix, weights=c, minlength=ny * nx)
        par_trame += np.bincount(tr, minlength=MAX_TRAMES)[:MAX_TRAMES]
        b = tr // fpb
        ok = b < nb
        blocs += np.bincount(b[ok] * (by * bx) + (y[ok] // BIN_DERIVE) * bx + x[ok] // BIN_DERIVE,
                             minlength=nb * by * bx)
        declin += np.bincount(b[ok] * config.CANAUX_TEMPS + c[ok], minlength=nb * config.CANAUX_TEMPS)
    trames_t = np.concatenate(trames_t)
    n_tr = len(trames_t)
    print(f"    {info['serie']}: {inten.sum():,} photons dans l'image, {n_tr} trames, "
          f"{nb} blocs de {fpb}  ({time.time() - t0:.0f} s)")
    return {"geo": geo, "nb": nb, "intensite": inten.reshape(ny, nx),
            "somme_canal": somme_c.reshape(ny, nx), "par_trame": par_trame[:n_tr],
            "trames_t": trames_t * geo.tic_s, "blocs": blocs.reshape(nb, by, bx),
            "declin_bloc": declin.reshape(nb, config.CANAUX_TEMPS)}


# ------------------------------------------------------------------ stimulation
def recaler_commandes(commandes, t_trames, par_trame):
    """Décalage δ (s) tel que t_cartes = t_NI + δ : maximise la corrélation entre la
    commande en créneau et la fuite de la stim dans ch2, image par image."""
    if not commandes:
        return 0.0, [], np.nan
    tc = t_trames + np.median(np.diff(t_trames)) / 2
    v = par_trame.astype(float)
    hp = v - ndi.median_filter(v, size=401, mode="nearest")     # retire la dérive lente
    hp /= np.std(hp) or 1.0
    best = (-np.inf, 0.0)
    for d in np.arange(-RECALAGE_MAX_S, RECALAGE_MAX_S + 1e-4, 0.005):
        ind = np.zeros_like(tc)
        for a, b, _ in commandes:
            ind[(tc >= a + d) & (tc < b + d)] = 1
        if ind.std() == 0:
            continue
        r = np.corrcoef(ind, hp)[0, 1]
        if r > best[0]:
            best = (r, d)
    if best[0] < MIN_CORR_FUITE:      # fuite trop faible (stim de basse puissance) : pas de recalage
        return 0.0, [(a, b) for a, b, _ in commandes], float(best[0])
    return float(best[1]), [(a + best[1], b + best[1]) for a, b, _ in commandes], float(best[0])


def fuite(fenetres, t_trames, par_trame):
    """Excédent de photons par image pendant chaque stim, en % de la base des 10 s avant."""
    out = []
    for a, b in fenetres:
        av = par_trame[(t_trames >= a - 10) & (t_trames < a)]
        pe = par_trame[(t_trames >= a + 0.5) & (t_trames < b - 0.5)]
        out.append(float((pe.mean() / av.mean() - 1) * 100) if len(av) and len(pe) else np.nan)
    return out


# ------------------------------------------------------------------ recalage
def vers_flim(im, ech, dy, dx, shape, order=1, cval=np.nan):
    """Image Bliq -> repère FLIM : x_flim = ech · x_bliq + dx, y_flim = y_bliq + dy."""
    return ndi.affine_transform(im.astype(float), [1, 1 / ech], offset=[-dy, -dx / ech],
                                order=order, output_shape=shape, cval=cval)


def recaler_fov(ref, im):
    """Meilleure (échelle x, dy, dx) de `im` (Bliq) sur `ref` (FLIM), et son NCC."""
    P = _prep(ref)
    m = np.zeros(ref.shape, bool)
    m[60:-60, 60:-60] = True
    best = None
    for e in ECHELLES_X:
        W = vers_flim(im, e, 0, 0, ref.shape, cval=0.0)
        Q = _prep(W)
        dy, dx = correlation_phase(P, Q)
        c = ncc(P, ndi.shift(Q, (dy, dx), order=0), m)
        if best is None or c > best[0]:
            best = (c, float(e), dy, dx)
    return {"ncc": best[0], "echelle_x": best[1], "dy": best[2], "dx": best[3]}


def derive(blocs, stim_bloc):
    """Décalage de chaque bloc par rapport à la moyenne des blocs hors stim, en pixels
    pleine résolution (corrélation de phase sous-pixel sur les images réduites)."""
    ref = _prep(blocs[~stim_bloc].mean(0), fin=1, large=12)
    m = np.zeros(ref.shape, bool)
    m[8:-8, 8:-8] = True
    out = []
    for b in blocs:
        Q = _prep(b, fin=1, large=12)
        dy, dx = correlation_phase(ref, Q, sous_pixel=True)
        out.append((dy * BIN_DERIVE, dx * BIN_DERIVE, ncc(ref, Q, m)))
    return np.array(out)


# ------------------------------------------------------------------ figures
def figure_par_image(dec, fen, out):
    """Champ entier image par image ; la dernière image, coupée par l'arrêt, est écartée."""
    t = dec[1]["trames_t"]
    n = len(t) - 1
    tm = (t + np.median(np.diff(t)) / 2)[:n]
    figures.figure_par_image(tm, dec[1]["par_trame"][:n], dec[2]["par_trame"][:n], fen,
                             "champ entier", out)


def figure_images(dec, fovs_alignes, modeles, out):
    s = dec[1]["intensite"] + dec[2]["intensite"]
    red = lambda im: im.reshape(256, 4, 256, 4).sum((1, 3))
    r = red(dec[2]["intensite"]) / np.maximum(red(dec[1]["intensite"]), 1)
    panneaux = [("FLIM ch1 (CFP), photons", dec[1]["intensite"], "gray"),
                ("FLIM ch2 (YFP), photons", dec[2]["intensite"], "gray"),
                ("ratio ch2 / ch1 (réduit 4 x 4)", r, "viridis")]
    for k in CANAUX:
        n = red(dec[k]["intensite"]).astype(float)
        tm = red(dec[k]["somme_canal"]) / np.maximum(n, 1) * modeles[k].dt - modeles[k].t_irf
        tm[n < 200] = np.nan
        panneaux.append((f"temps moyen ch{k} − IRF (ns, 4 x 4)", tm, "magma"))
    for nom, im in fovs_alignes.items():
        panneaux.append((f"Bliq {nom} (C2, recalée)", im, "gray"))
    nc = 4
    nl = int(np.ceil(len(panneaux) / nc))
    fig, axes = plt.subplots(nl, nc, figsize=(4.6 * nc, 4.8 * nl), facecolor=SURFACE)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (titre, im, cmap) in zip(axes.ravel(), panneaux):
        h = ax.imshow(im, cmap=cmap, interpolation="nearest",
                      **contraste(im, 2 if cmap != "gray" else 1, 98 if cmap != "gray" else 99.7))
        ax.set_title(titre, color=INK, fontsize=10)
        if cmap != "gray":
            cb = fig.colorbar(h, ax=ax, fraction=0.046, pad=0.02)
            cb.ax.tick_params(labelsize=7, colors=INK_2, length=0)
    fig.suptitle(f"Images de la session — {s.sum():,.0f} photons, ch1 + ch2", color=INK,
                 fontsize=12.5, x=0.01, ha="left")
    fig.tight_layout()
    sauver(fig, out, dpi=110)


def figure_derive(t_bloc, d, fen, out):
    fig, ax = plt.subplots(figsize=(10, 3.6), facecolor=SURFACE)
    style(ax)
    bandes_stim(ax, fen)
    ax.plot(t_bloc, d[:, 0], color="#2a78d6", lw=1.5, marker="o", ms=2.5, label="dy")
    ax.plot(t_bloc, d[:, 1], color="#eb6834", lw=1.5, marker="o", ms=2.5, label="dx")
    ax.axhline(0, color=INK_2, lw=0.8, ls="--")
    ax.set_ylabel("décalage (px)", color=INK)
    ax.set_xlabel("temps (s)", color=INK)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(f"Dérive du champ, bloc par bloc (max |d| = {np.abs(d[:, :2]).max():.2f} px ; "
                 f"NCC {d[:, 2].min():.2f}–{d[:, 2].max():.2f})", color=INK, fontsize=11, loc="left")
    fig.tight_layout()
    sauver(fig, out)


def figure_recalage(ref, fovs_alignes, tr, out):
    n = len(fovs_alignes)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.4), facecolor=SURFACE, squeeze=False)
    norm = lambda im: np.clip((im - np.nanpercentile(im, 1)) /
                              (np.nanpercentile(im, 99.5) - np.nanpercentile(im, 1)), 0, 1)
    a = norm(ref)
    for ax, (nom, im) in zip(axes[0], fovs_alignes.items()):
        b = norm(np.nan_to_num(im, nan=0.0))
        ax.imshow(np.dstack([a, b, a]), interpolation="nearest")
        t = tr[nom]
        plan = "même plan" if t["meme_plan"] else "AUTRE PLAN FOCAL"
        ax.set_title(f"{nom} ({t['date'][11:16]}) — {plan}\nNCC avec le FLIM {t['ncc_flim']:.2f}, "
                     f"avec {t['reference']} {t['ncc_reference']:.2f}", color=INK, fontsize=9.5)
        ax.axis("off")
    t = next(iter(tr.values()))
    fig.suptitle(f"Recalage Bliq -> FLIM (échelle x {t['echelle_x']:.4f}, dy {t['dy']:+.0f}, "
                 f"dx {t['dx']:+.0f}, commun à toutes les images) : magenta = FLIM (ch1 + ch2), "
                 f"vert = Bliq C2, blanc = superposition", color=INK, fontsize=11.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    sauver(fig, out, dpi=110)


# ------------------------------------------------------------------ main
def main():
    root = config.racine()
    D = config.dirs(root)
    S = config.session(root)
    fpb = S["images_par_bloc"]
    print(f"{root.name} — session du {S['debut']}, blocs de {fpb} images")

    dec = {k: decoder_canal(S["spc"][k], fpb) for k in CANAUX}
    nb = min(dec[k]["nb"] for k in CANAUX)
    t_tr = dec[1]["trames_t"]
    t_bloc0 = t_tr[np.arange(nb) * fpb]
    t_bloc1 = t_tr[np.minimum((np.arange(nb) + 1) * fpb, len(t_tr) - 1)]

    # --- fenêtres de stim, recalées sur l'horloge des cartes
    delta, fen, r = recaler_commandes(S["commandes"], t_tr, dec[2]["par_trame"][:len(t_tr)])
    g = config.STIM_GARDE_S
    stim_bloc = np.zeros(nb, bool)
    for a, b in fen:                    # toutes les fenêtres masquent leurs blocs…
        stim_bloc |= (t_bloc1 > a - g) & (t_bloc0 < b + g)
    # … mais une stim coupée par la fin de l'enregistrement (ou bien plus courte que les
    # autres) n'a pas d'« après » mesurable : elle n'entre pas dans l'analyse.
    if fen:
        d_med = np.median([b - a for a, b in fen])
        coupees = [(a, b) for a, b in fen if b - a < 0.5 * d_med or b > t_tr[-1] - 0.5]
        fen = [w for w in fen if w not in coupees]
    else:
        coupees = []
    fuites = {k: fuite(fen, t_tr + np.median(np.diff(t_tr)) / 2, dec[k]["par_trame"][:len(t_tr)])
              for k in CANAUX}
    print(f"  stimulation : {len(fen)} fenêtres, décalage NI -> cartes {delta * 1000:+.0f} ms "
          f"(corrélation {r:.2f}"
          + (", trop faible : horloge NI gardée)" if r < MIN_CORR_FUITE else ")"))
    for a, b in coupees:
        print(f"    stim incomplète {a:7.2f}–{b:7.2f} s : masquée, mais écartée de l'analyse")
    for i, (a, b) in enumerate(fen):
        print(f"    stim {i + 1} : {a:7.2f}–{b:7.2f} s   fuite ch1 {fuites[1][i]:+.1f} %, "
              f"ch2 {fuites[2][i]:+.1f} %")

    # --- contrôle : temps des blocs contre frames.csv
    fr = config.frames_flimapp(root)
    p0 = {int(x["pass"]): float(x["pass_start_s"]) for x in fr}
    ecart = [abs(p0[b + 1] - t_bloc0[b]) for b in range(nb) if b + 1 in p0]
    print(f"  début des blocs contre frames.csv : écart max {max(ecart) * 1000:.3f} ms "
          f"sur {len(ecart)} blocs" if ecart else "  frames.csv vide")

    # --- dérive
    blocs = (dec[1]["blocs"][:nb] + dec[2]["blocs"][:nb]).astype(float)
    d = derive(blocs, stim_bloc)
    print(f"  dérive : max |dy| {np.abs(d[:, 0]).max():.2f} px, max |dx| {np.abs(d[:, 1]).max():.2f} px")

    # --- images fixes Bliq
    # Le passage Bliq -> FLIM est instrumental (même scanner, deux horloges de pixel) :
    # on le mesure sur l'image qui ressemble le plus au FLIM et on l'applique à toutes.
    # Chaque image est ensuite comparée à cette référence : en dessous de
    # MIN_NCC_MEME_PLAN, elle n'est pas au plan focal de la session.
    ref = (dec[1]["intensite"] + dec[2]["intensite"]).astype(float)
    tr, alignes, brutes = {}, {}, {}
    for f in config.fovs(root):
        if "C2" not in f["canaux"]:
            continue
        nom = f"{f['nm']} nm {f['moment']}"
        im = tifffile.imread(f["canaux"]["C2"]).squeeze().astype(float)
        t = recaler_fov(ref, im)
        tr[nom] = {"dossier": f["nom"], "nm": f["nm"], "moment": f["moment"],
                   "date": config.ome_date(f["canaux"]["C2"]), "ncc_propre": t["ncc"],
                   "propre": {k: t[k] for k in ("echelle_x", "dy", "dx")}}
        brutes[nom] = (f, im)
    if tr:
        nref = max(tr, key=lambda n: tr[n]["ncc_propre"])
        commun = tr[nref]["propre"]
        P = _prep(ref)
        m = np.zeros(ref.shape, bool)
        m[60:-60, 60:-60] = True
        R = _prep(vers_flim(brutes[nref][1], commun["echelle_x"], commun["dy"], commun["dx"],
                            ref.shape, cval=0.0))
        print(f"  passage Bliq -> FLIM mesuré sur {nref} : échelle x {commun['echelle_x']:.4f}, "
              f"dy {commun['dy']:+.0f}, dx {commun['dx']:+.0f}")
        for nom, (f, im) in brutes.items():
            Q = _prep(vers_flim(im, commun["echelle_x"], commun["dy"], commun["dx"], ref.shape, cval=0.0))
            t = tr[nom]
            t.update(commun, reference=nref, ncc_flim=ncc(P, Q, m), ncc_reference=ncc(R, Q, m))
            t["meme_plan"] = bool(nom == nref or t["ncc_reference"] >= MIN_NCC_MEME_PLAN)
            print(f"  {f['nom']:24s} {t['date'][11:16]}  NCC FLIM {t['ncc_flim']:.2f}  "
                  f"NCC {nref} {t['ncc_reference']:.2f}  "
                  f"{'même plan' if t['meme_plan'] else 'AUTRE PLAN FOCAL'}")
            for c, p in f["canaux"].items():
                im_c = tifffile.imread(p).squeeze().astype(float)
                alignes[f"{nom} {c}"] = vers_flim(im_c, commun["echelle_x"], commun["dy"],
                                                  commun["dx"], ref.shape)
    order = sorted(tr, key=lambda n: (tr[n]["moment"] != "avant", tr[n]["nm"]))
    tr = {n: tr[n] for n in order}

    # --- sorties
    I, dt = config.irf(root)
    modeles = {k: Modele(I[k], dt) for k in CANAUX}
    np.savez_compressed(
        D["projections"] / "flim.npz",
        **{f"intensite_ch{k}": dec[k]["intensite"].astype(np.uint32) for k in CANAUX},
        **{f"somme_canal_ch{k}": dec[k]["somme_canal"].astype(np.float32) for k in CANAUX},
        **{f"par_trame_ch{k}": dec[k]["par_trame"] for k in CANAUX},
        **{f"blocs_ch{k}": dec[k]["blocs"][:nb].astype(np.uint16) for k in CANAUX},
        **{f"declin_bloc_ch{k}": dec[k]["declin_bloc"][:nb] for k in CANAUX},
        trames_t=t_tr, somme=ref.astype(np.float32))
    np.savez_compressed(D["projections"] / "fovs.npz",
                        **{k.replace(" ", "_"): v.astype(np.float32) for k, v in alignes.items()})
    info = {
        "jeu": root.name, "debut": S["debut"], "images_par_bloc": fpb, "n_blocs": int(nb),
        "periode_image_s": float(np.median(np.diff(t_tr))),
        "blocs": [{"bloc": b + 1, "t0": round(float(t_bloc0[b]), 4), "t1": round(float(t_bloc1[b]), 4),
                   "stim": bool(stim_bloc[b])} for b in range(nb)],
        "commandes_NI": [{"t0": a, "t1": b, "valeur": v} for a, b, v in S["commandes"]],
        "decalage_NI_cartes_s": delta, "correlation_fuite": r,
        "stims": [{"t0": round(a, 3), "t1": round(b, 3), "fuite_ch1_pct": fuites[1][i],
                   "fuite_ch2_pct": fuites[2][i]} for i, (a, b) in enumerate(fen)],
        "stims_incompletes": [{"t0": round(a, 3), "t1": round(b, 3)} for a, b in coupees],
        "derive": [{"bloc": b + 1, "dy": float(d[b, 0]), "dx": float(d[b, 1]), "ncc": float(d[b, 2])}
                   for b in range(nb)],
        "fovs": tr,
        "photons": {f"ch{k}": int(dec[k]["intensite"].sum()) for k in CANAUX},
    }
    (D["projections"] / "session.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))

    figure_par_image(dec, fen, D["qc"] / "par_image.png")
    figure_derive((t_bloc0 + t_bloc1) / 2, d, fen, D["qc"] / "derive.png")
    sel = {n.replace(" C2", ""): v for n, v in alignes.items() if n.endswith("C2")}
    if sel:
        figure_recalage(ref, sel, tr, D["qc"] / "recalage_fovs.png")
    figure_images(dec, sel, modeles, D["projections"] / "images.png")
    print(f"\nSorties -> {D['projections']}")


if __name__ == "__main__":
    main()
