"""Lance les étapes du pipeline sur un ou plusieurs jeux de données.

Usage :  python3 executer.py <jeu> [<jeu> ...] [--de 1] [--a 7] [--diameter 75] [--image racine]
                             [--masques-de <jeu du même champ>]
ex.      python3 executer.py ../SCLM+eNpHR_2026-10-09/*/ --de 2 --diameter 75

- l'étape 2 (Cellpose) tourne avec ~/.venvs/cellpose/bin/python, les autres avec le
  python3 courant ; --diameter, --image et --masques-de lui sont passés ;
- l'étape 3 n'écrit la sélection par défaut (tous les ROIs gardés) que s'il n'y en a pas
  encore : un tri manuel n'est jamais écrasé ;
- un jeu sans ligne de base (1re stim en tout début de session) s'arrête proprement aux
  étapes 5 à 7 (valeurs brutes seulement).
"""
import argparse
import subprocess
import sys
from pathlib import Path

ICI = Path(__file__).resolve().parent
CELLPOSE_PY = Path.home() / ".venvs" / "cellpose" / "bin" / "python"
ETAPES = {1: "step1_decode.py", 2: "step2_cellpose.py", 3: "step3_select.py",
          4: "step4_declins.py", 5: "step5_analyse.py", 6: "step6_ratio_corrige.py",
          7: "step7_retour.py"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jeux", nargs="+")
    ap.add_argument("--de", type=int, default=1)
    ap.add_argument("--a", type=int, default=7)
    ap.add_argument("--diameter", type=float)
    ap.add_argument("--image", choices=("racine", "brut"))
    ap.add_argument("--masques-de", dest="source")
    args = ap.parse_args()
    for jeu in map(Path, args.jeux):
        if not (jeu / "data" / "FLIM").exists():
            print(f"!! {jeu} : pas de data/FLIM, ignoré")
            continue
        print(f"\n######## {jeu.resolve().name}")
        for n in range(args.de, args.a + 1):
            cmd = [str(CELLPOSE_PY if n == 2 else sys.executable), str(ICI / ETAPES[n]), str(jeu)]
            if n == 2:
                cmd += (["--diameter", str(args.diameter)] if args.diameter else []) + \
                       (["--image", args.image] if args.image else []) + \
                       (["--masques-de", args.source] if args.source else [])
            if n == 3:
                if (jeu / "selection" / "cellpose_selection.json").exists():
                    print("  étape 3 : sélection existante gardée")
                    continue
                cmd.append("--tout")
            print(f"  --- étape {n}")
            r = subprocess.run(cmd, cwd=ICI)
            if r.returncode and n < 5:
                sys.exit(f"!! étape {n} en échec sur {jeu}")


if __name__ == "__main__":
    main()
