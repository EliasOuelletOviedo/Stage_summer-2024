"""Étape 1 — recalage des conditions sur Basal, puis image somme par FOV.

Pour chaque FOV : les 2 canaux (C2, C3) des 3 conditions sont recalés sur Basal
par corrélation de phase, puis additionnés en une seule image (6 images -> 1).
Sortie : PNG + TIFF de la somme, et un .npz des canaux recalés pour l'étape 4.
"""
import glob
import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import DATA_ROOT, CONDITIONS, CHANNELS, FOV_MAP, SLICES, N_FOV, fov_key, dirs

MAX_SHIFT = 60      # dérive de platine plausible entre conditions (~30 µm)
MIN_NCC = 0.15      # sous ce score le pic de corrélation n'est pas fiable


def load_channel(slice_, cond, fov, ch):
    hits = sorted(glob.glob(str(DATA_ROOT / slice_ / cond / fov / "Bliq VMS" / ch / "*.tif")))
    if not hits:
        raise FileNotFoundError(f"{slice_}/{cond}/{fov}/{ch}")
    return tifffile.imread(hits[0]).squeeze().astype(np.float64)


def _prep(im):
    """Passe-bande + z-score : rend comparables des images d'intensités très différentes."""
    im = gaussian_filter(im, 1.5) - gaussian_filter(im, 20.0)
    return (im - im.mean()) / (im.std() + 1e-9)


def estimate_shift(moving, reference):
    """Translation entière (dy, dx) qui amène `moving` sur `reference`, + score NCC."""
    a, b = _prep(reference), _prep(moving)
    cc = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real / a.size)
    c = np.array(cc.shape) // 2
    w = cc[c[0] - MAX_SHIFT:c[0] + MAX_SHIFT + 1, c[1] - MAX_SHIFT:c[1] + MAX_SHIFT + 1]
    i = np.unravel_index(np.argmax(w), w.shape)
    return (int(i[0] - MAX_SHIFT), int(i[1] - MAX_SHIFT)), float(w.max())


def apply_shift(im, dy, dx, fill=np.nan):
    """Translation entière sans interpolation ; les bords entrants valent `fill`."""
    out = np.full_like(im, fill, dtype=np.float64)
    ys, xs = im.shape
    y0d, y1d = max(dy, 0), min(ys + dy, ys)
    y0s, y1s = max(-dy, 0), min(ys - dy, ys)
    x0d, x1d = max(dx, 0), min(xs + dx, xs)
    x0s, x1s = max(-dx, 0), min(xs - dx, xs)
    out[y0d:y1d, x0d:x1d] = im[y0s:y1s, x0s:x1s]
    return out


def process_fov(slice_, idx, D, log):
    key = fov_key(slice_, idx)
    ref = None
    reg, valid = {}, {}

    for cond in CONDITIONS:
        fov = FOV_MAP[slice_][cond][idx]
        if fov is None:
            log.append(f"{key}\t{cond}\t--\t\t\t\tABSENT")
            continue
        chans = {ch: load_channel(slice_, cond, fov, ch) for ch in CHANNELS}
        composite = sum(chans.values())

        if cond == "Basal":
            ref, (dy, dx), score, note = composite, (0, 0), 1.0, "ref"
        else:
            (dy, dx), score = estimate_shift(composite, ref)
            if score < MIN_NCC:
                # Image trop bruitée : un pic de corrélation faible donne un
                # décalage arbitraire, bien pire que pas de recalage du tout.
                dy, dx, note = 0, 0, "NON RECALE (NCC faible)"
            else:
                note = "ok"

        for ch in CHANNELS:
            reg[f"{cond}_{ch}"] = apply_shift(chans[ch], dy, dx)
        v = apply_shift(np.ones_like(composite), dy, dx, fill=0.0) > 0.5
        valid[cond] = v
        log.append(f"{key}\t{cond}\t{fov}\tdy={dy:+d}\tdx={dx:+d}\tNCC={score:.3f}\t{note}")

    # Somme des 2 canaux x 3 conditions, bords non recouverts comptés comme 0.
    total = np.zeros((1024, 1024))
    for arr in reg.values():
        total += np.nan_to_num(arr, nan=0.0)

    np.savez_compressed(D["registered"] / f"{key}.npz",
                        **{k: v.astype(np.float32) for k, v in reg.items()},
                        **{f"valid_{c}": m for c, m in valid.items()},
                        summed=total.astype(np.float32))
    tifffile.imwrite(D["sum_images"] / f"{key}_sum.tif", total.astype(np.float32))

    lo, hi = np.percentile(total, (1, 99.7))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(total, cmap="gray", vmin=lo, vmax=hi)
    ax.set_title(f"{key} — somme C2+C3 x (Basal, High_K, Return)", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(D["sum_images"] / f"{key}_sum.png", dpi=130)
    plt.close(fig)
    return key


def main():
    D = dirs()
    log = ["fov\tcondition\tdossier\tdy\tdx\tNCC\tstatut"]
    keys = []
    for slice_ in SLICES:
        for idx in range(N_FOV):
            keys.append(process_fov(slice_, idx, D, log))
            print(f"  ok  {keys[-1]}")
    (D["qc"] / "registration_log.tsv").write_text("\n".join(log) + "\n")
    print(f"\n{len(keys)} images somme -> {D['sum_images']}")


if __name__ == "__main__":
    main()
