"""Recalage par corrélation de phase (translation, sous-pixel) et ressemblance (NCC)."""
import numpy as np
from scipy import ndimage as ndi


def _prep(im, sigma=1.0):
    a = ndi.gaussian_filter(np.asarray(im, float), sigma)
    a = (a - a.mean()) / (a.std() + 1e-12)
    w = np.hanning(a.shape[0])[:, None] * np.hanning(a.shape[1])[None, :]
    return a * w


def decalage(im, ref, sigma=1.0):
    """(dy, dx) tel que im(r) ≈ ref(r − d) : le contenu de `im` est déplacé de d par rapport
    à `ref` (les masques de `ref` se déplacent de +d pour suivre `im`). Sous-pixel par
    parabole autour du pic."""
    A, B = np.fft.fft2(_prep(im, sigma)), np.fft.fft2(_prep(ref, sigma))
    R = A * np.conj(B)
    r = np.real(np.fft.ifft2(R / (np.abs(R) + 1e-12)))
    i = np.unravel_index(np.argmax(r), r.shape)
    d = []
    for ax, n in enumerate(r.shape):
        idx = list(i)
        c = r[tuple(idx)]
        idx[ax] = (i[ax] - 1) % n
        m = r[tuple(idx)]
        idx[ax] = (i[ax] + 1) % n
        p = r[tuple(idx)]
        den = m - 2 * c + p
        sub = 0.5 * (m - p) / den if den != 0 else 0.0
        v = i[ax] + sub
        d.append(v if v < n / 2 else v - n)
    return float(d[0]), float(d[1])


def ncc(a, b, masque=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if masque is not None:
        a, b = a[masque], b[masque]
    a, b = a - a.mean(), b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-300))


def deplacer(im, d, ordre=0):
    """Déplace le contenu de `im` de d = (dy, dx) pixels entiers (bords : 0)."""
    dy, dx = int(round(d[0])), int(round(d[1]))
    return ndi.shift(im, (dy, dx), order=ordre, cval=0)
