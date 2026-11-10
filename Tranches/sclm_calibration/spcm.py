"""Lecture des flux .spc écrits par le logiciel SPCM avec la SPC-QC-104 (mode FIFO_IMAGE).

Un en-tête de 4 octets (bits 0-23 : macrotemps en fs, 2 048 131 fs = 2,048131 ns), puis des
enregistrements de 32 bits petit-boutistes (format établi par FLIMApp, `format_fifo_qc104.jl`,
vérifié contre la DLL de B&H) :

    bit 31 seul (& 0xC0000000 == 0x80000000)  un débordement du macrotemps (4096 tics)
    bit 30                                     marqueur : M0-M3 sur les bits 12-15
                                               (M1 = début de ligne, M2 = début de trame)
    bits 31-30 nuls                            photon : macrotemps 0-11, routage 12-15,
                                               canal du TDC 16-27 (4 ps, temps croissant
                                               depuis la synchro), entrée 28-29 (0 = IN1)

Rangement des photons : une ligne commence à son M1, la première ligne qui suit un M2 est la
ligne 0 de la trame ; ligne de l'image = rang dans la trame − DECALAGE_LIGNES, colonne =
⌊(t − t_M1) · tic / dwell⌋ − DECALAGE_PIXELS (le pixel de Bliq, 51,82 ns). Le microtemps
(canal de 4 ps) est tiré uniformément dans sa largeur, ramené modulo la période du laser
puis rangé en 256 canaux de P/256 : sans ce tirage, 12 ou 13 canaux du TDC tombent dans
chaque canal de 48,7 ps (peigne de ±4 %), comme dans FLIMApp et pileup/.
"""
import numpy as np

import config as C

ENTETE_FS = 2048131


def lire(path, graine=0):
    """Photons d'un cycle, rangés dans l'image.

    Rend un dict : y, x (int16), entree (int8 : 0 = IN1 = YFP, 1 = IN2 = CFP), canal
    (uint8, 0..255 sur la période), t_s (float32, s depuis le premier événement) pour les
    photons de l'image ; passages (512,) : nombre de balayages complets de chaque ligne de
    l'image ; n_trames ; duree_s ; n_entree (photons de chaque entrée, image ou non).
    """
    raw = np.fromfile(path, dtype="<u4")
    if len(raw) < 2:
        raise ValueError(f"{path} : fichier vide")
    if int(raw[0]) & 0xFFFFFF != ENTETE_FS:
        raise ValueError(f"{path} : macrotemps inattendu ({int(raw[0]) & 0xFFFFFF} fs)")
    w = raw[1:]
    del raw
    deb = (w & 0xC0000000) == 0x80000000
    ov = np.cumsum(deb, dtype=np.int64)            # débordements jusqu'à chaque mot
    del deb
    marq = (w & 0x40000000) != 0
    ph = (w & 0xC0000000) == 0

    def temps(sel):
        return ov[sel] * 4096 + (w[sel] & 0xFFF).astype(np.int64)

    hm = (w[marq] >> 12) & 0xF
    tm = temps(marq)
    L = tm[(hm & 2) != 0]                          # M1 : débuts de ligne
    F = tm[(hm & 4) != 0]                          # M2 : débuts de trame
    tp = temps(ph)
    wp = w[ph]
    del ov, w, marq, ph
    t0 = min(tp[0] if len(tp) else np.inf, tm[0] if len(tm) else np.inf)
    t1 = max(tp[-1] if len(tp) else -np.inf, tm[-1] if len(tm) else -np.inf)
    entree = ((wp >> 28) & 3).astype(np.int8)
    n_entree = np.bincount(entree, minlength=4)[:2]

    # ligne et trame de chaque photon
    trL = np.searchsorted(F, L, side="right") - 1          # trame de chaque ligne (-1 : avant le 1er M2)
    premier = np.searchsorted(trL, trL, side="left")
    rang = np.arange(len(L)) - premier                     # rang de la ligne dans sa trame
    k = np.searchsorted(L, tp, side="right") - 1           # dernière ligne commencée avant le photon
    kk = np.maximum(k, 0)
    y = rang[kk] - C.DECALAGE_LIGNES
    x = np.floor((tp - L[kk]) * (C.TIC_NS / C.DWELL_NS)).astype(np.int64) - C.DECALAGE_PIXELS
    ok = (k >= 0) & (trL[kk] >= 0) & (y >= 0) & (y < C.LIGNES_IMAGE) & (x >= 0) & (x < C.PIXELS_LIGNE)

    # passages : lignes complètes (suivies d'un autre M1) de chaque rang de l'image
    complete = np.r_[np.ones(max(len(L) - 1, 0), bool), np.zeros(min(len(L), 1), bool)] & (trL >= 0)
    yr = rang[complete] - C.DECALAGE_LIGNES
    yr = yr[(yr >= 0) & (yr < C.LIGNES_IMAGE)]
    passages = np.bincount(yr, minlength=C.LIGNES_IMAGE).astype(np.int32)

    # microtemps -> canal de la période
    rng = np.random.default_rng(graine)
    tdc = ((wp[ok] >> 16) & 0xFFF).astype(np.float64)
    micro = np.mod((tdc + rng.random(len(tdc))) * C.CANAL_TDC_NS, C.PERIODE_NS)
    canal = np.minimum((micro / C.PERIODE_NS * C.N_CANAUX).astype(np.int32), C.N_CANAUX - 1).astype(np.uint8)
    return dict(y=y[ok].astype(np.int16), x=x[ok].astype(np.int16), entree=entree[ok], canal=canal,
                t_s=((tp[ok] - t0) * C.TIC_NS * 1e-9).astype(np.float32),
                passages=passages, n_trames=int(len(F)), duree_s=float((t1 - t0) * C.TIC_NS * 1e-9),
                n_entree=n_entree.astype(np.int64))


def image(ph):
    """(2, 512, 1024) photons par pixel, canal 0 = CFP, 1 = YFP (ordre de config.CANAUX)."""
    out = np.zeros((2, C.LIGNES_IMAGE, C.PIXELS_LIGNE), np.uint16)
    for i, nom in enumerate(C.CANAUX):
        s = ph["entree"] == C.ENTREE[nom]
        out[i] = np.bincount(ph["y"][s].astype(np.int64) * C.PIXELS_LIGNE + ph["x"][s],
                             minlength=C.LIGNES_IMAGE * C.PIXELS_LIGNE).reshape(C.LIGNES_IMAGE, C.PIXELS_LIGNE)
    return out


def declin_total(ph, nom):
    """Déclin (256 canaux) de tous les photons de l'image d'un canal."""
    return np.bincount(ph["canal"][ph["entree"] == C.ENTREE[nom]], minlength=C.N_CANAUX)


def declins_rois(ph, masques, n_rois, bins_s=None):
    """Déclins par ROI : (n_rois + 1, 2, 256) — ligne 0 = hors des ROIs — et, si `bins_s`
    est donné (nombre de tranches de 1 s), photons par ROI et par seconde (n_rois + 1, 2, bins_s)."""
    lab = masques[ph["y"], ph["x"]].astype(np.int64)
    dec = np.zeros((n_rois + 1, 2, C.N_CANAUX), np.int64)
    sec = np.zeros((n_rois + 1, 2, bins_s), np.int64) if bins_s else None
    for i, nom in enumerate(C.CANAUX):
        s = ph["entree"] == C.ENTREE[nom]
        dec[:, i] = np.bincount(lab[s] * C.N_CANAUX + ph["canal"][s],
                                minlength=(n_rois + 1) * C.N_CANAUX).reshape(n_rois + 1, C.N_CANAUX)
        if bins_s:
            b = np.minimum(ph["t_s"][s].astype(np.int64), bins_s - 1)
            sec[:, i] = np.bincount(lab[s] * bins_s + b, minlength=(n_rois + 1) * bins_s).reshape(n_rois + 1, bins_s)
    return dec, sec
