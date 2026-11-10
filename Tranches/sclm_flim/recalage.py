"""Recalage d'images par corrélation de phase — commun aux étapes 1 (dérive, images
fixes Bliq) et 2 (ROIs repris d'une autre session du même champ)."""
import numpy as np
from scipy import ndimage as ndi


def _prep(im, fin=2, large=30):
    im = ndi.gaussian_filter(im.astype(float), fin)
    return im - ndi.gaussian_filter(im, large)


def correlation_phase(A, B, sous_pixel=False):
    """Décalage (dy, dx) à appliquer à B pour qu'il recouvre A."""
    R = np.fft.fft2(A) * np.conj(np.fft.fft2(B))
    R /= np.abs(R) + 1e-12
    c = np.real(np.fft.ifft2(R))
    i = np.unravel_index(np.argmax(c), c.shape)
    out = []
    for ax, (v, n) in enumerate(zip(i, c.shape)):
        d = float(v if v < n // 2 else v - n)
        if sous_pixel:
            idx = [list(i), list(i)]
            idx[0][ax] = (v - 1) % n
            idx[1][ax] = (v + 1) % n
            cm, cp, c0 = c[tuple(idx[0])], c[tuple(idx[1])], c[i]
            den = cm - 2 * c0 + cp
            if den < 0:
                d += 0.5 * (cm - cp) / den
        out.append(d)
    return tuple(out)


def ncc(A, B, m):
    a, b = A[m] - A[m].mean(), B[m] - B[m].mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-300))
