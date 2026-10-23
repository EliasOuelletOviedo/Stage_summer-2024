"""Lecture des fichiers B&H : flux .spc des deux cartes, .sdt (déclins, IRF, taux).

Flux .spc — un en-tête de 4 octets, puis des mots de 32 bits petit-boutistes.

SPC-150N (FIFO_150) : bits 0-11 macrotemps, 12-15 routage ou marqueurs, 16-27 ADC,
28 MARK, 29 GAP, 30 MTOV, 31 INVALID. INVALID + MTOV sans MARK : plusieurs débordements
(nombre dans les bits 0-27) ; MARK + INVALID : marqueur ; INVALID seul : photon détecté
hors des limites du TAC (il compte pour le temps mort, pas pour les déclins). Le
macrotemps (25 ns ici) est dans les bits 0-23 de l'en-tête, en dixièmes de ns. Le temps
croît quand l'ADC décroît : canal = 4095 − ADC.

SPC-QC-104 (format de FLIMApp, `format_fifo_qc104.jl`) : bits 0-11 macrotemps (2,048 ns),
12-15 routage ou marqueurs, 16-27 canal du TDC (4 ps, temps croissant), 28-29 entrée ;
mot & 0xC0000000 == 0x80000000 : un débordement (4096 tics) ; bit 30 : marqueur ;
bits 30-31 nuls : photon.

Marqueurs : M1 = début de ligne, M2 = début de trame.
"""
from dataclasses import dataclass

import numpy as np
import sdtfile

from config import CANAL_QC_NS, N_CANAUX, PERIODE_NS, PIXEL_NS, PIXELS_LIGNE, TIC_QC_S


@dataclass
class Flux:
    """Photons d'une entrée : temps absolu (ns), canal de la période (0..255, −1 = sans
    microtemps), et marqueurs de ligne / trame (ns)."""
    t_ns: np.ndarray
    canal: np.ndarray
    lignes_ns: np.ndarray
    trames_ns: np.ndarray
    micro_ns: np.ndarray        # microtemps fin (ns) ; NaN si inconnu
    resolution_t_ns: float      # pas du macrotemps


def _marqueurs(t, est_marq, haut):
    L = t[est_marq & (((haut >> 1) & 1) == 1)]
    F = t[est_marq & (((haut >> 2) & 1) == 1)]
    return L, F


def lire_150n(path):
    raw = np.fromfile(path, dtype="<u4")
    tic_ns = (int(raw[0]) & 0xFFFFFF) * 0.1
    w = raw[1:]
    inv = (w >> 31) & 1
    mtov = (w >> 30) & 1
    mark = (w >> 28) & 1
    multi = (inv == 1) & (mtov == 1) & (mark == 0)
    inc = np.where(multi, (w & 0x0FFFFFFF).astype(np.int64) * 4096, mtov.astype(np.int64) * 4096)
    t = np.cumsum(inc) + (w & 0xFFF).astype(np.int64)
    haut = ((w >> 12) & 0xF).astype(np.int8)
    est_marq = ~multi & (mark == 1) & (inv == 1)
    photon = ~multi & (mark == 0)                    # valide ou hors limites du TAC
    valide = photon & (inv == 0)
    L, F = _marqueurs(t, est_marq, haut)
    adc = ((w >> 16) & 0xFFF).astype(np.int32)
    canal = np.where(valide, (4095 - adc) >> 4, -1)[photon]
    micro = np.where(valide, (4095 - adc + 0.5) * 12.5 / 4096, np.nan)[photon]
    return Flux(t[photon] * tic_ns, canal.astype(np.int16), L * tic_ns, F * tic_ns,
                micro, tic_ns)


def lire_qc104(path, graine=0):
    """{entrée (1, 2): Flux}. Le canal de la période vient du microtemps du TDC, ramené
    sur [0, P) (0,01 % des photons dépassent à peine P). Comme dans FLIMApp, le canal de
    4 ps est tiré uniformément dans sa largeur avant d'être rangé en canaux de P/256 :
    sans ce tirage, 12 ou 13 canaux du TDC tombent dans chacun (peigne de ±4 %)."""
    w = np.fromfile(path, dtype="<u4")[1:]
    tic_ns = TIC_QC_S * 1e9
    deb = (w & 0xC0000000) == 0x80000000
    marq = (w & 0x40000000) != 0
    photon = (w & 0xC0000000) == 0
    t = np.cumsum(deb.astype(np.int64) * 4096) + (w & 0xFFF).astype(np.int64)
    haut = ((w >> 12) & 0xF).astype(np.int8)
    L, F = _marqueurs(t, marq, haut)
    entree = ((w >> 28) & 3)[photon]
    u = np.random.default_rng(graine).random(int(photon.sum()))
    micro = (((w >> 16) & 0xFFF)[photon] + u) * CANAL_QC_NS
    micro = np.mod(micro, PERIODE_NS)
    canal = np.minimum((micro / PERIODE_NS * N_CANAUX).astype(np.int16), N_CANAUX - 1)
    tp = t[photon] * tic_ns
    return {e + 1: Flux(tp[entree == e], canal[entree == e], L * tic_ns, F * tic_ns,
                        micro[entree == e], tic_ns) for e in (0, 1)}


def lire_spc(carte, dossier, puissance):
    """{canal: Flux} d'une acquisition."""
    if carte == "SPC-QC-104":
        return lire_qc104(dossier / f"{puissance}.spc")
    return {k: lire_150n(dossier / f"{puissance}_m{k}.spc") for k in (1, 2)}


# ------------------------------------------------------------------ image
def pixels(f: Flux):
    """(y, x, exposition) : position de chaque photon dans l'image (−1 hors image) et nombre
    de passages du balayage sur chaque ligne y. Ligne = rang depuis le M2 de sa trame,
    colonne = temps depuis le M1 / temps par pixel ; pas de décalage : seule la carte
    d'intensité sert, pour le taux local."""
    L, F = f.lignes_ns, f.trames_ns
    k = np.searchsorted(L, f.t_ns, side="right") - 1
    tr_L = np.searchsorted(F, L, side="right") - 1            # trame de chaque ligne
    premier = np.searchsorted(tr_L, tr_L, side="left")
    rang_L = np.arange(len(L)) - premier
    ok_L = tr_L >= 0
    kk = np.maximum(k, 0)
    x = np.floor((f.t_ns - L[kk]) / PIXEL_NS).astype(np.int64)
    y = rang_L[kk]
    # La dernière ligne n'est complète que si un M1 la suit : ses photons sont écartés.
    ok = (k >= 0) & (k < len(L) - 1) & ok_L[kk] & (x >= 0) & (x < PIXELS_LIGNE)
    complete = ok_L[:-1]
    ny = int(rang_L[:-1][complete].max()) + 1 if complete.any() else 1
    expo = np.bincount(rang_L[:-1][complete], minlength=ny).astype(float)
    return np.where(ok, y, -1), np.where(ok, x, -1), expo


# ------------------------------------------------------------------ .sdt
def lire_sdt(path):
    return sdtfile.SdtFile(path)


def info_mesure(s, i=0):
    """Taux de StopInfo (moyennes sur la mesure) et durée."""
    m = s.measure_info[i]
    m = m.flatten()[0] if hasattr(m, "flatten") else m
    si = m["StopInfo"]
    return {"duree_s": float(si["stop_time"]), "sync": float(si["max_sync_rate"]),
            "cfd": float(si["max_cfd_rate"]), "tac": float(si["max_tac_rate"]),
            "adc": float(si["max_adc_rate"]), "carte": m["mod_type"].decode()}


def irf(dossier, carte, canal):
    """IRF sur la grille des déclins (256 canaux sur une période), fond retiré, somme 1.

    QC-104 : histogramme de 64 ps sur 16,384 ns, rééchantillonné sur P/256 par sa
    fonction cumulée. 150N : déjà 256 canaux sur 12,5 ns."""
    from config import IRF
    nom, bloc, ligne = IRF[carte][canal]
    s = lire_sdt(dossier / nom)
    y = s.data[bloc][ligne].astype(float)
    t = s.times[bloc] * 1e9
    dt = t[1] - t[0]
    pic = int(np.argmax(y))
    loin = np.abs(np.arange(len(y)) - pic) > 40
    if carte == "SPC-QC-104":
        loin &= t < PERIODE_NS - dt                          # au-delà de P : pas de photons
    y = np.clip(y - np.median(y[loin]), 0, None)
    y[np.abs(np.arange(len(y)) - pic) > 40] = 0              # le fond n'est pas l'IRF
    if carte == "SPC-QC-104":
        bords = np.r_[t, t[-1] + dt]
        cum = np.r_[0, np.cumsum(y)]
        g = np.arange(N_CANAUX + 1) * PERIODE_NS / N_CANAUX
        y = np.diff(np.interp(g, bords, cum))
    return y / y.sum()
