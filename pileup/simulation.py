"""Monte-Carlo de l'empilement : des photons poissoniens par impulsion, un temps mort.

Chaque impulsion laser envoie au détecteur un nombre de photons de Poisson, de moyenne
μ_det ; leur temps dans la période suit le `gabarit` (le déclin mesuré au plus bas taux,
étendu aux canaux où la carte ne mesure pas). Un photon n'est enregistré que si la carte
n'est pas morte : temps mort non paralysable, qui commence au photon enregistré.

- QC-104 : temps mort `d` par entrée, ~7,4 ns, mesuré sur les intervalles entre photons
  (étape 1). Il est plus court que la période : un 2e photon de la même impulsion peut
  être enregistré s'il arrive plus de d après le 1er ; et un photon tardif rend la carte
  aveugle au début de la période suivante.
- 150N : temps mort ~100 ns, plusieurs périodes ; deux variantes :
  « photon » : D compté depuis le photon (la carte se libère à la phase φ + D mod P) ;
  « synchro » : la conversion part de la synchro qui suit le photon (TAC en start-stop
  inversé) et dure D : la carte se libère toujours à la même phase.

Le premier photon d'une impulsion est privilégié (empilement : déclin plus rapide) ;
une carte qui se libère au milieu d'une période ne voit que la fin de cette période
(déclin plus lent). Le résultat net dépend du temps mort et de la forme du déclin.
"""
import numpy as np

from config import N_CANAUX, PERIODE_NS


def gabarit(declin, fenetre):
    """Déclin étendu à toute la période : hors de la fenêtre de la carte, le niveau moyen
    des 5 canaux voisins (fond, queue de l'impulsion précédente). Somme 1."""
    g = np.asarray(declin, float).copy()
    lo, hi = fenetre
    g[:lo] = g[lo:lo + 5].mean()
    g[hi:] = g[hi - 5:hi].mean()
    return g / g.sum()


def simuler(gab, mu_det, temps_mort_ns, mode="photon", photons=1_500_000, graine=0):
    """Histogramme enregistré (256 canaux), μ enregistré et μ détecté par impulsion.

    `photons` : nombre de photons arrivant au détecteur à simuler (taille de l'échantillon) ;
    le nombre d'impulsions en découle. `mode` : « photon » ou « synchro » (voir en tête)."""
    rng = np.random.default_rng(graine)
    n_imp = max(int(photons / mu_det), 1)
    n = rng.poisson(mu_det * n_imp)
    k = np.sort(rng.integers(0, n_imp, n))                 # impulsion de chaque photon
    canal = rng.choice(N_CANAUX, size=n, p=gab)
    phase = (canal + rng.random(n)) * PERIODE_NS / N_CANAUX
    T = k * PERIODE_NS + phase
    ordre = np.argsort(T, kind="stable")
    T, canal, k = T[ordre], canal[ordre], k[ordre]
    gardes = []
    i = 0
    while i < n:
        gardes.append(i)
        fin = T[i] + temps_mort_ns if mode == "photon" else (k[i] + 1) * PERIODE_NS + temps_mort_ns
        i = int(np.searchsorted(T, fin, side="left"))
        if i <= gardes[-1]:
            i = gardes[-1] + 1
    gardes = np.asarray(gardes)
    h = np.bincount(canal[gardes], minlength=N_CANAUX).astype(float)
    return h, len(gardes) / n_imp, n / n_imp


def flux_poisson(gab, mu_det, photons=1_500_000, graine=0):
    """(temps absolus en ns, canal, nombre d'impulsions) des photons arrivant au détecteur,
    triés : le même flux sert à comparer plusieurs temps morts (courbes lisses en D)."""
    rng = np.random.default_rng(graine)
    n_imp = max(int(photons / mu_det), 1)
    n = rng.poisson(mu_det * n_imp)
    k = np.sort(rng.integers(0, n_imp, n))
    canal = rng.choice(N_CANAUX, size=n, p=gab)
    T = k * PERIODE_NS + (canal + rng.random(n)) * PERIODE_NS / N_CANAUX
    ordre = np.argsort(T, kind="stable")
    return T[ordre], canal[ordre], n_imp


def garder(T, temps_mort_ns, T_liste=None):
    """Photons enregistrés avec un temps mort non paralysable compté depuis chaque photon
    enregistré (masque). Un photon à plus de D du précédent est toujours gardé : seuls ceux
    qui en sont plus proches sont examinés un à un. `T_liste` : T.tolist(), à passer quand
    le même flux sert à plusieurs D (la conversion coûte plus que la boucle)."""
    Tl = T.tolist() if T_liste is None else T_liste
    garde = np.ones(len(Tl), bool)
    proches = (np.flatnonzero(np.diff(T) < temps_mort_ns) + 1).tolist()
    jetes = []
    dernier, precedent = 0.0, -2
    for i in proches:
        if i - 1 != precedent:            # le photon d'avant était loin du sien : gardé
            dernier = Tl[i - 1]
        if Tl[i] - dernier < temps_mort_ns:
            jetes.append(i)
        else:
            dernier = Tl[i]
        precedent = i
    garde[jetes] = False
    return garde


def temps_mort_effectif(mu_enr, mu_det):
    """Temps mort de la formule non paralysable, 1/r_enr − 1/r_det, en ns."""
    return (1 / mu_enr - 1 / mu_det) * PERIODE_NS
