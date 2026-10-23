"""Configuration partagée pour l'analyse SCLM 2026-08-17 (tranches, ratio C2/C3)."""
from pathlib import Path

DATA_ROOT = Path("/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-08-17_tranches/SCLM")
OUT_ROOT = Path(__file__).resolve().parent.parent / "SCLM_2026-08-17"

CONDITIONS = ["Basal", "High_K", "Return"]
REF_CONDITION = "Basal"          # toutes les conditions sont recalées sur Basal
CHANNELS = ["C2", "C3"]          # C1 ignoré
PIXEL_SIZE_UM = 0.5087830142943353
CELL_DIAMETER_PX = 30            # soma ~15 µm à 0.509 µm/px

# Appariement des FOVs, vérifié par corrélation croisée invariante en translation.
# Slice 3 : Return n'a que 4 acquisitions ; Return-004 correspond à Basal-005,
# et Basal-004 n'a aucun Return (None).
FOV_MAP = {
    "Slice 1": {
        "Basal":  [f"FOV_{i}-001" for i in range(1, 6)],
        "High_K": [f"FOV_{i}-001" for i in range(1, 6)],
        "Return": [f"FOV-00{i}" for i in range(1, 6)],
    },
    "Slice 2": {
        "Basal":  [f"FOV-00{i}" for i in range(1, 6)],
        "High_K": [f"FOV-00{i}" for i in range(1, 6)],   # FOV-006 = timelapse, exclu
        "Return": [f"FOV-00{i}" for i in range(1, 6)],
    },
    "Slice 3": {
        "Basal":  [f"FOV-00{i}" for i in range(1, 6)],
        "High_K": [f"FOV-00{i}" for i in range(1, 6)],
        "Return": ["FOV-001", "FOV-002", "FOV-003", None, "FOV-004"],
    },
}
SLICES = list(FOV_MAP)
N_FOV = 5


def fov_key(slice_, idx):
    """Identifiant court et stable d'un FOV, ex. 'Slice1_FOV4'."""
    return f"{slice_.replace(' ', '')}_FOV{idx + 1}"


def dirs():
    d = {k: OUT_ROOT / k for k in
         ("sum_images", "registered", "masks", "selection", "results", "qc")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d
