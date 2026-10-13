"""Configuration partagée — SCLM + eNpHR, tranches 2026-09-21, ratio FRET C2/C3 dans le temps.

Même logique que `../sclm/` : ROIs triés à la main, puis ratio C2/C3 en valeurs brutes.
La différence est que chaque acquisition est ici un **décours temporel** (C2 et C3
entrelacés, T impairs = C2, T pairs = C3) et non une image fixe.
"""
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "Data"

# Dossier de sortie. SCLM_OUT permet d'écrire ailleurs sans toucher aux résultats du
# pipeline principal (c'est ce que fait train_FOV11.py).
_out = os.environ.get("SCLM_OUT")
OUT_ROOT = (Path(_out) if Path(_out).is_absolute() else ROOT / _out) if _out else ROOT

# FOVs analysés à part, par leur propre script : le pipeline principal les ignore.
# SCLM_FOVS (liste séparée par des virgules) restreint l'analyse à ces FOVs-là.
SEPARATE_FOVS = ("FOV11",)

CHANNELS = ["C2", "C3"]              # ratio FRET = C2/C3, comme dans ../sclm
PIXEL_SIZE_UM = 0.5087830142943353   # métadonnée OME
CELL_DIAMETER_PX = 60                # soma ~30 µm à 0.509 µm/px (mesuré sur les sommes)

MIN_FRAMES_TIMELAPSE = 5   # en dessous : image fixe (800/950/1064 nm), pas de décours

# Images fixes : seul C2 est analysé, et seulement à 800 et 950 nm. C3 y est quasiment
# vide hors 800 nm, et 1064 nm n'est pas retenu.
FIXED_CHANNEL = "C3"
FIXED_WAVELENGTHS = (800, 950)

# --- détection de la fenêtre de stimulation ------------------------------------
# La lumière de stim fuit dans C2 sur tout le champ (de +11 % à x25 le niveau de base
# selon la puissance) : ces images ne mesurent plus le FRET.
#
# La détection porte sur les **bords** du plateau dans la moyenne pleine image de C2,
# pas sur son niveau : ce qui distingue une stim d'un photoblanchiment, c'est que
# l'allumage et l'extinction tiennent en une image. Un seuil sur le niveau échouait dès
# que le plateau occupait une bonne part de l'acquisition, puisqu'il gonflait lui-même
# l'estimation du bruit (jusqu'à x25 sur FOV3) et se hissait hors de portée.
#
# Le bord est mesuré en comparant les moyennes de STIM_EDGE_W images de part et
# d'autre, et non deux images voisines : une marche franche répond de toute sa
# hauteur, une dérive lente de pente x W seulement, et le bruit est divisé par
# sqrt(W/2). C'est ce qui rend détectable FOV5/stim_VIS_1-001, dont la fuite n'est que
# de +2,9 % (3,3 écarts sur des différences brutes, 9,1 avec la fenêtre).
STIM_EDGE_W = 5        # demi-fenêtre du contraste de bord, en images
STIM_STEP_K = 6.0      # bord minimal, en écarts robustes de ce contraste
STIM_MIN_RISE = 0.01   # ... et plateau au moins 1 % au-dessus de la base
STIM_MIN_FRAMES = 10   # ... et durant au moins 10 images. Doit rester > STIM_EDGE_W :
                       # la fenêtre étale une image aberrante isolée sur W images, ce
                       # qui imiterait un court plateau.
STIM_PAD = 1           # images de garde de chaque côté : l'ouverture et la fermeture de
                       # l'obturateur tombent en milieu d'image, qui est alors à moitié éclairée

# --- recalage ------------------------------------------------------------------
BLOCK = 25        # images par bloc pour la dérive intra-acquisition (une seule
                  # image est trop bruitée pour un pic de corrélation fiable)
MAX_SHIFT = 60    # dérive de platine plausible (~30 µm)
MIN_NCC = 0.15    # sous ce score le pic n'est pas fiable : aucun décalage appliqué

MIN_VALID_FRAC = 0.90   # part du ROI qui doit rester dans le champ après recalage

# Fenêtre de réponse immédiate, en secondes après la fin de la stim. La moyenne sur
# toute la période d'après est sensible à la dérive lente de la ligne de base ; les
# premières secondes le sont bien moins.
RESPONSE_WINDOW_S = 10.0
PRE_WINDOW_S = 10.0     # fenêtre avant la stim pour la figure avant / après

# Étape 5 : les N_BEST_TRACES meilleures réponses (ROI x stim), toutes données
# confondues. Une acquisition n'est candidate que si elle laisse au moins MIN_POST_S
# après la stim : écarte FOV1/stim-002 (88 images, 29 s après la stim), dont la durée et
# le calendrier diffèrent des autres.
N_BEST_TRACES = 50
MIN_POST_S = 45.0


def _t_index(p):
    """Numéro T du nom de fichier, ex. 'stim-001-C2-T137.tif' -> 137."""
    return int(re.search(r"-T(\d+)\.tif$", p.name).group(1))


def slug(name):
    """Nom de fichier sûr : '800 nm-001' -> '800nm-001'."""
    return re.sub(r"[^A-Za-z0-9_.-]", "", name.replace(" ", ""))


def frame_files(fov, rec, ch):
    """Images d'un canal d'une acquisition, triées par numéro T croissant."""
    d = DATA_ROOT / fov / rec / "Bliq VMS" / ch
    return sorted(d.glob("*.tif"), key=_t_index)


def wavelength(rec):
    """Longueur d'onde d'une image fixe, lue en tête de nom : '800 nm-001' -> 800."""
    m = re.match(r"\s*(\d+)", rec)
    return int(m.group(1)) if m else None


def fixed_recordings(recs):
    """Images fixes retenues (C2 à 800 / 950 nm), par longueur d'onde croissante.

    `recs` est une liste (acquisition, n_images, horodatage) telle que la donne survey().
    """
    keep = [r for r, n, _ in recs
            if n < MIN_FRAMES_TIMELAPSE and wavelength(r) in FIXED_WAVELENGTHS]
    return sorted(keep, key=wavelength)


def stim_type(rec):
    """Type de stimulation d'après le nom : 'VIS', 'IR', ou None si non étiqueté
    (FOV1/stim-001, stim-002 et FOV2/stim_max-001)."""
    if "IR" in rec:
        return "IR"
    if "VIS" in rec:
        return "VIS"
    return None


def acquisition_date(fov, rec):
    """Horodatage OME de la première image — sert à ordonner les acquisitions."""
    import tifffile
    with tifffile.TiffFile(frame_files(fov, rec, "C2")[0]) as tf:
        xml = tf.pages[0].tags[270].value
    return re.search(r"<AcquisitionDate>(.*?)</", xml).group(1)


def delta_t(path):
    """DeltaT (s depuis le début de l'acquisition) de l'en-tête OME, sans décoder les pixels."""
    import tifffile
    with tifffile.TiffFile(path) as tf:
        xml = tf.pages[0].tags[270].value
    return float(re.search(r'DeltaT="([\d.eE+-]+)"', xml).group(1))


def survey():
    """Contenu de Data/ : {FOV: [(acquisition, n_images, horodatage), ...]} en ordre chronologique.

    Une acquisition n'est retenue que si C2 et C3 existent et ont le même nombre d'images.
    """
    only = os.environ.get("SCLM_FOVS")
    only = set(only.split(",")) if only else None
    out = {}
    for fov_dir in sorted(p for p in DATA_ROOT.iterdir() if p.is_dir()):
        if (fov_dir.name not in only) if only else (fov_dir.name in SEPARATE_FOVS):
            continue
        recs = []
        for rec_dir in sorted(p for p in fov_dir.iterdir() if p.is_dir()):
            rec = rec_dir.name
            n = [len(frame_files(fov_dir.name, rec, ch)) for ch in CHANNELS]
            if not all(n) or n[0] != n[1]:
                continue
            recs.append((rec, n[0], acquisition_date(fov_dir.name, rec)))
        if recs:
            out[fov_dir.name] = sorted(recs, key=lambda r: r[2])
    return out


def reference_recording(recs):
    """Acquisition servant de repère spatial au FOV : le plus long décours, puis le plus ancien."""
    return min(recs, key=lambda r: (-r[1], r[2]))[0]


def dirs():
    d = {k: OUT_ROOT / k for k in ("projections", "masks", "selection", "results", "qc")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d
