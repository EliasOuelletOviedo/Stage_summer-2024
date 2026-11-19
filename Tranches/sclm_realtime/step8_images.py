"""Étape 8 — reconstruire les images sans les marqueurs de ligne et de trame (M1/M2).

Usage :  python3 step8_images.py <dossier>          (~2 min ; après l'étape 1)

En Realtime par passes, FLIMApp n'enregistre pas M1/M2 : un photon n'a que son temps. Le
polygone balaie pourtant le champ (ligne de 55,564 µs, trame de 540 lignes) et le débit de
photons porte la structure de l'image. La phase des lignes, elle, dérive par rapport à
l'horloge de la carte (~20 µs/s, avec des ondulations) : on la suit trame par trame.

  1. période nominale T0 de la ligne (raie à 17 997 Hz du débit de photons) ; fenêtres d'une
     trame (540·T0 ≈ 30 ms) ; dans chacune, le profil de la ligne : photons repliés modulo T0
     sur NB canaux (tous les photons des deux canaux, phases opto et pauses comprises) — sur
     une trame entière, c'est la projection de l'image sur x, la même d'une trame à l'autre ;
  2. décalage de chaque trame par corrélation croisée circulaire avec un gabarit (sous-canal
     par interpolation parabolique), le gabarit refait des profils réalignés (3 itérations) ;
     décalages déroulés d'une trame à l'autre (une ligne glissée compte : la trame en dépend) ;
  3. chaque photon : u = (t − t₀)/T0·NB − s(t), s interpolé linéairement entre les centres
     des trames ; ligne L = ⌊u / NB⌋, x = u − L·NB, y = L mod 540. L'origine de y est
     arbitraire : l'image est roulée pour mettre le trait brillant du retour du
     galvanomètre (la ligne la plus saillante sur son voisinage) en bas ;
  4. recalage : une image par bloc de BLOC_S secondes, recalée sur l'image de toute la session
     (corrélation croisée 2D, sous-pixel) ; ce décalage corrige à la fois une dérive de la
     tranche et une erreur lente du suivi (elles ne se distinguent pas). Image recalée = somme
     des blocs ramenés sur la référence.

Contrôles : écart-type des décalages entre trames consécutives ; décalages des blocs au fil de
la session ; pas de motif de facette (décalage x ligne à ligne, vérifié le 2026-11-11 :
0,26 canal).

Sorties : mesures/images.npz, results/images.png, qc/images_suivi.png
"""
import time

import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.ndimage import shift as decaler

import config as C
import figures as G
import spc as S
from figures import plt

F_LIGNE_HZ = 17997.345     # raie du débit de photons (2026-11-11), affinée par le suivi
NB = 1072                  # canaux par ligne : ~51,8 ns, le pixel du Bliq
LIGNES = 540               # lignes par trame (réglage du scanner)
ITERATIONS = 3
BLOC_S = 10.0              # blocs recalés sur l'image de la session


def decalage(a, b):
    """s (canaux, sous-canal) tel que a ≈ roll(b, s), par corrélation croisée circulaire."""
    c = np.real(np.fft.ifft(np.fft.fft(a) * np.conj(np.fft.fft(b))))
    k = int(np.argmax(c))
    y0, y1, y2 = c[k - 1], c[k], c[(k + 1) % len(c)]
    den = y0 - 2 * y1 + y2
    s = k + (0.5 * (y0 - y2) / den if den != 0 else 0.0)
    return s - len(a) if s > len(a) / 2 else s


def recalage(a, b):
    """(dy, dx) sous-pixel tel que a ≈ b décalé de (dy, dx), par corrélation croisée 2D."""
    c = np.real(np.fft.ifft2(np.fft.fft2(a - a.mean()) * np.conj(np.fft.fft2(b - b.mean()))))
    ky, kx = np.unravel_index(np.argmax(c), c.shape)
    out = []
    for k, n, ax in ((ky, c.shape[0], 0), (kx, c.shape[1], 1)):
        y0 = c[(ky - 1) % c.shape[0], kx] if ax == 0 else c[ky, (kx - 1) % c.shape[1]]
        y2 = c[(ky + 1) % c.shape[0], kx] if ax == 0 else c[ky, (kx + 1) % c.shape[1]]
        y1 = c[ky, kx]
        den = y0 - 2 * y1 + y2
        v = k + (0.5 * (y0 - y2) / den if den != 0 else 0.0)
        out.append(v - n if v > n / 2 else v)
    return tuple(out)


def flux(ses):
    """(temps en s, canal 0 = CFP / 1 = YFP, code, canal de temps, est_photon) par bloc, les deux fichiers fusionnés
    bloc par bloc n'est pas nécessaire : on lit un fichier après l'autre."""
    for k, nom in enumerate(("cfp", "yfp")):
        for t, ph, _, champ, adc, _ in S._blocs(ses["spc"][C.CANAL[nom]]["spc"]):
            yield k, t[ph] * ses["tic_s"], champ[ph], (4095 - adc[ph]) // 16


def main():
    root = C.racine()
    D = C.dirs(root)
    ses = C.session(root)
    P = C.charger_passes(D)
    t0s = time.time()
    T0 = 1 / F_LIGNE_HZ
    W = LIGNES * T0
    A0 = float(P["t_debut_s"][0]) - 1.0
    A1 = float(P["t_fin_s"][-1]) + 1.0
    nf = int((A1 - A0) / W) + 1

    # --- 1. profils par trame
    prof = np.zeros((nf, NB))
    for _, t, _, _ in flux(ses):
        m = (t >= A0) & (t < A1)
        t = t[m]
        j = ((t - A0) / W).astype(np.int64)
        x = (((t - A0) / T0) % 1.0 * NB).astype(np.int64)
        prof += np.bincount(j * NB + x, minlength=nf * NB).reshape(nf, NB)
    n_par_trame = prof.sum(1)
    garde = n_par_trame > 0.25 * np.median(n_par_trame)
    print(f"  {nf} trames de {1e3 * W:.2f} ms ; {int(n_par_trame.sum())} photons ({time.time() - t0s:.0f} s)")

    # --- 2. décalages
    tpl = prof[garde][:20].sum(0)
    for it in range(ITERATIONS):
        sh = np.array([decalage(prof[j], tpl) if garde[j] else np.nan for j in range(nf)])
        ok = np.isfinite(sh)
        d = np.diff(sh[ok])
        d = (d + NB / 2) % NB - NB / 2
        shu = np.full(nf, np.nan)
        shu[ok] = np.r_[sh[ok][0], sh[ok][0] + np.cumsum(d)]
        tpl = np.zeros(NB)
        for j in np.flatnonzero(ok):
            tpl += np.roll(prof[j], -int(round(sh[j])))
        dd = np.diff(shu[ok])
        print(f"  itération {it + 1} : écart-type des différences entre trames {dd.std():.2f} canaux "
              f"({1e9 * dd.std() * T0 / NB:.0f} ns), saut max {np.abs(dd).max():.1f}")
    tc = A0 + (np.arange(nf) + 0.5) * W
    tc, shu = tc[ok], shu[ok]
    derive = (shu[-1] - shu[0]) / NB / (tc[-1] - tc[0])          # lignes par seconde
    f_vraie = F_LIGNE_HZ - derive

    # --- 3. images
    def pixels(t):
        u = (t - A0) / T0 * NB - np.interp(t, tc, shu)
        L = np.floor(u / NB).astype(np.int64)
        return L, np.minimum((u - L * NB).astype(np.int64), NB - 1)

    L0 = None
    nbl = int(np.ceil((tc[-1] - tc[0]) / BLOC_S))
    img = np.zeros((2, 2, LIGNES, NB))           # [canal CFP/YFP, tous / acquisition seule]
    blk = np.zeros((nbl, LIGNES, NB), np.float32)  # blocs de BLOC_S s, tous les photons
    for k, t, code, _ in flux(ses):
        m = (t >= tc[0]) & (t < tc[-1])
        t, code = t[m], code[m]
        L, x = pixels(t)
        if L0 is None:
            L0 = int(L.min())
        y = (L - L0) % LIGNES
        idx = y * NB + x
        img[k, 0] += np.bincount(idx, minlength=LIGNES * NB).reshape(LIGNES, NB)
        acq = code != 0
        img[k, 1] += np.bincount(idx[acq], minlength=LIGNES * NB).reshape(LIGNES, NB)
        b = ((t - tc[0]) / BLOC_S).astype(np.int64)
        for bb in np.unique(b):
            blk[bb] += np.bincount(idx[b == bb], minlength=LIGNES * NB).reshape(LIGNES, NB)
    # roulement : le trait brillant du retour du galvanomètre (ligne la plus saillante) en bas
    rang = img[:, 0].sum(0).sum(1)
    k = int(np.argmax(rang - gaussian_filter1d(rang, 6, mode="wrap")))
    roul = -(k + 4)                                # la ligne k + 4 devient la première
    img = np.roll(img, roul, axis=2)
    blk = np.roll(blk, roul, axis=1)
    # recalage des blocs sur l'image de la session
    ref = gaussian_filter(img[:, 0].sum(0), 2)
    dep = np.array([recalage(gaussian_filter(b, 2), ref) for b in blk])          # (dy, dx) par bloc
    t_bloc = tc[0] + (np.arange(nbl) + 0.5) * BLOC_S
    rec = np.zeros_like(ref)
    for b, (dy, dx) in zip(blk, dep):
        rec += decaler(b, (-dy, -dx), order=1, mode="wrap")
    print(f"  fréquence de ligne suivie : {f_vraie:.4f} Hz ; dérive de phase moyenne {1e6 * derive * T0:.1f} µs/s ; "
          f"recalage des blocs de {BLOC_S:.0f} s : dy de {dep[:, 0].min():+.1f} à {dep[:, 0].max():+.1f}, "
          f"dx de {dep[:, 1].min():+.1f} à {dep[:, 1].max():+.1f} px")
    # images recalées par canal : chaque bloc ramené sur la référence, au prorata des photons du canal
    np.savez_compressed(D["mesures"] / "images.npz", img=img.astype(np.float32), recale=rec.astype(np.float32),
                        t_bloc=t_bloc, dep_bloc=dep, bloc_s=BLOC_S, tc=tc, decalages=shu, T0=T0, NB=NB, A0=A0, L0=L0,
                        roulement=roul, lignes=LIGNES, gabarit=tpl, f_ligne=f_vraie)
    figures(D, img, rec, tc, shu, tpl, T0, t_bloc, dep)
    print(f"  -> mesures/images.npz ({time.time() - t0s:.0f} s)")


def figures(D, img, rec, tc, shu, tpl, T0, t_bloc, dep):
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    tot = img[:, 0].sum(0)
    for ax, (im, nom) in zip(axes.ravel(), ((rec, "tous les photons (CFP + YFP), toute la session, blocs recalés"),
                                           (img[1, 0], "YFP"), (img[0, 0], "CFP"),
                                           (gaussian_filter(img[1, 0], 3) / np.maximum(gaussian_filter(img[0, 0], 3), 1),
                                            "rapport YFP / CFP (lissé 3 px)"))):
        if nom.startswith("rapport"):
            msk = gaussian_filter(tot, 3) > np.percentile(gaussian_filter(tot, 3), 60)
            v = np.where(msk, im, np.nan)
            h = ax.imshow(v, cmap="viridis", vmin=np.nanpercentile(v, 2), vmax=np.nanpercentile(v, 98))
            fig.colorbar(h, ax=ax, fraction=0.025)
        else:
            g = gaussian_filter(im, 1.0)
            ax.imshow(g, cmap="gray", vmin=0, vmax=np.percentile(g, 99.7))
        ax.set_title(nom, color=G.INK, fontsize=11, loc="left")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Images reconstruites sans M1/M2 (phase des lignes suivie trame par trame) — 2026-11-11",
                 color=G.INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    G.sauver(fig, D["results"] / "images.png", dpi=110)

    fig, axes = plt.subplots(1, 3, figsize=(20, 4.8))
    for ax in axes:
        G.style(ax)
    axes[0].plot(tc, 1e6 * (shu - shu[0]) * T0 / NB, color=G.INK, lw=0.8)
    axes[0].set_xlabel("temps (s)", color=G.INK)
    axes[0].set_ylabel("décalage de phase des lignes (µs)", color=G.INK)
    G.titre(axes[0], "Phase du polygone par rapport à l'horloge de la carte", "une estimation par trame (30 ms)", fontsize=10.5)
    d = np.diff(shu) * T0 / NB * 1e9
    axes[1].hist(d - np.median(d), bins=120, color=G.TRACE)
    axes[1].set_xlabel("différence entre trames consécutives, médiane retirée (ns)", color=G.INK)
    G.titre(axes[1], "Bruit du suivi", f"écart-type {np.std(d):.0f} ns (pixel : 51,8 ns)", fontsize=10.5)
    axes[2].plot(t_bloc, dep[:, 1], "o-", color=G.INK, ms=3, lw=0.8, label="dx")
    axes[2].plot(t_bloc, dep[:, 0], "o-", color=G.PALETTE[1], ms=3, lw=0.8, label="dy")
    axes[2].set_xlabel("temps (s)", color=G.INK)
    axes[2].set_ylabel("décalage du bloc (px)", color=G.INK)
    axes[2].legend(frameon=False, fontsize=9)
    G.titre(axes[2], "Recalage des blocs sur l'image de la session", "dérive de la tranche et erreur lente du suivi",
            fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["qc"] / "images_suivi.png")


if __name__ == "__main__":
    main()
