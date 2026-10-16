#=
===============================================================================
 Simulateur de stimulation de tornades par faisceau laser
===============================================================================

 Variante deux lasers. A chaque position la tornade est balayee deux fois de
 suite sans pause : d'abord avec le laser 1, puis avec le laser 2. Un delai
 suit, les deux lasers eteints, avant de sauter a la position suivante.
 La sequence complete peut etre repetee. Le resultat est une video MP4.

 Modele de temps
 ---------------
 Chaque passe occupe DUREE_TORNADE_MS et compte N_POINTS points : le faisceau
 s'arrete sur chaque point pendant DUREE_TORNADE_MS / N_POINTS, exactement
 comme le buffer de points envoye aux galvos. Une position occupe donc
 2 x DUREE_TORNADE_MS + DELAI_MS. Pendant DELAI_MS les deux lasers sont
 eteints et la position saute instantanement au premier point de la suivante.

 Rendu
 -----
 1 pixel = 1 mV par defaut, donc la resolution suit la plage des galvos
 (+/- 1000 mV en x et en y donne 2001 x 2001 pixels). Le spot est une
 gaussienne d'ecart-type SIGMA_FAISCEAU. Toute la tornade est peinte d'un coup
 pendant sa passe, puis s'efface d'un coup : l'image montre la dose integree
 d'une passe entiere, pas la position instantanee du faisceau.

 Lancement : julia Simulateurs/Tornades.jl
===============================================================================
=#

using Colors
using Dates
using Random
using Printf
using ProgressBars
using FFMPEG_jll

# =============================================================================
#  PARAMETRES
# =============================================================================

# --- Geometrie d'une tornade -------------------------------------------------
const N_POINTS = 1000         # nombre de points par tornade ; la cadence de
                              # points reste celle d'avant, la passe etant
                              # huit fois plus courte
const N_TOURS  = 4.0          # nombre de tours complets de la spirale
                              # N_POINTS / N_TOURS donne les points par tour :
                              # sous environ 2*pi*RAYON/SIGMA_FAISCEAU points
                              # par tour, la spirale se lit en pointilles.
                              # INTERPOLER_TRAJET la rend continue.
const RAYON    = 200.0        # rayon exterieur de la tornade [mV]

# --- Sequence temporelle -----------------------------------------------------
# Cadencement cale sur le budget d'images. Une seconde de video vaut dix
# secondes d'experience, donc a 30 images par seconde une image couvre
# 10000/30 ms. Une passe en occupe deux, la pause une seule, donc une ROI vaut
# 5 images et une sequence de 4 ROI en vaut 20.
const DUREE_TORNADE_MS = 20000 / 30  # 2 images par passe laser
const DELAI_MS         = 10000 / 30  # 1 image de pause
const N_REPETITIONS    = 90          # 90 sequences x 20 images = 1800 images

# --- Positions des tornades --------------------------------------------------
# Le nombre de tornades est la longueur de cette liste. Elles sont visitees
# dans l'ordre ou elles sont ecrites.
const POSITIONS = [
    (-500.0,  500.0),
    ( 500.0,  500.0),
    (-500.0, -500.0),
    ( 500.0, -500.0),
]

# --- Faisceau ----------------------------------------------------------------
const SIGMA_FAISCEAU = 15.0                 # ecart-type du spot gaussien [mV]
const COULEUR_LASER_1 = colorant"#ff0000"
const COULEUR_LASER_2 = colorant"#ffb300"

# Ordre de passage. false : laser 1 puis laser 2. true : laser 2 puis laser 1.
# Les couleurs, les graphiques et les sens de spirale restent attaches a leur
# laser, seul l'ordre des deux passes change.
const INVERSER_ORDRE_LASERS = false

# Sens radial de la spirale, reglable pour chaque laser separement.
# true  : du centre vers le bord (sortante)
# false : du bord vers le centre (rentrante)
# Dans les deux cas la rotation reste anti-horaire.
const SPIRALE_SORTANTE_1 = true
const SPIRALE_SORTANTE_2 = false

# --- Champ de vue et resolution ----------------------------------------------
const XLIM = (-1000.0, 1000.0)   # plage des galvos en x [mV]
const YLIM = (-1000.0, 1000.0)   # plage des galvos en y [mV]
const ECHELLE_PX = 1.0           # pixels par mV ; 1.0 donne 2001 x 2001 pour
                                 # +/- 1000 mV. Baisser a 0.5 pour un apercu
                                 # deux fois plus rapide.

# --- Video -------------------------------------------------------------------
const FPS     = 30.0            # images par seconde de la video
const RALENTI = 0.1             # 1 s de video = 10 s d'experience
const DOSSIER_SORTIE = joinpath(@__DIR__, "sorties")
const NOM_FICHIER    = ""        # vide = nom horodate automatique

# --- Concentration, asservissement et clampage -------------------------------
# Chaque tornade du laser 1 constitue une acquisition : elle rend une valeur de
# [Cl] pour sa ROI. Le laser 2 est l'actionneur, sa puissance fait monter la
# concentration. Entre deux acquisitions la commande est bloquee, et le systeme
# du premier ordre est integre exactement sur le pas.
#
#   tau_i d[Cl]_i/dt = -([Cl]_i - CL_BASAL[i]) + GAIN_COMMANDE * u_i
#
# Une valeur par ROI, dans le meme ordre que POSITIONS.
const CL_BASAL  = [8.0, 16.0, 12.0, 14.0]          # concentration basale [mM]
# Le retour vers la basale n'est pas asservi : il ne depend que de tau. Entrer
# dans la bande des 2 % demande ln(1/0.02) = 3.91 tau, d'ou un retour de 211 a
# 282 s selon la ROI, 250 s en moyenne, soit 25 s de video.
const TAU_ROI_S = [54.0, 62.0, 68.0, 72.0]       # constante de temps du premier ordre [s]

const CL_CIBLE      = 28.0   # consigne commune a toutes les ROI [mM]
const GAIN_COMMANDE = 30.0   # elevation de [Cl] en regime etabli a u = 1 [mM]

# Gains choisis par balayage pour ces CL_BASAL, TAU_ROI_S et CL_CIBLE. Ils
# donnent une montee deux fois plus lente que l'optimum en temps, afin que la
# forme exponentielle du premier ordre soit lisible : etablissement en 14 s au
# lieu de 7, et ecart de 1.6 % seulement a une exponentielle pure.
# Pour la reponse la plus rapide possible, prendre Kp = 0.32 et Ki = 0.14 : la
# commande reste alors collee a la butee sur toute la montee, qui se lit comme
# une rampe et non comme une exponentielle.
# A refaire si tau, le gain de commande ou la cadence d'acquisition changent.
# Gains choisis par balayage sous trois contraintes : la commande ne touche
# jamais la butee, aucun depassement, et 95 % de l'echelon atteints entre 133 et
# 153 s selon la ROI, soit 13 a 15 s de video. La montee reste ainsi une
# exponentielle franche.
# Le correcteur tourne une fois par acquisition et l'integrale accumule e*dt :
# si toutes les durees sont multipliees par dix, il faut diviser Ki par dix pour
# retrouver exactement la meme suite discrete. Kp, lui, ne bouge pas.
const PI_KP = 0.0445         # gain proportionnel [commande / mM]
const PI_KI = 0.00068        # gain integral [commande / (mM s)]

# Le clampage est cale sur les numeros d'acquisition, la cadence du correcteur.
const CLAMP_DEBUT_ACQ = 16   # montee a 100 s (periode d'acquisition 20/3 s)
const CLAMP_FIN_ACQ   = 46   # retour a 300 s, exponentiel vers la basale

const BRUIT_MESURE = 0.0     # ecart-type du bruit de mesure [mM], 0 = aucun
const GRAINE_BRUIT = 1234    # graine, pour que deux rendus soient identiques

const CL_PLOT = (0.0, 36.0)  # bornes verticales du graphe de concentration

# Palette categorielle des ROI, une couleur par serie, jamais recyclee. Valeurs
# validees pour la vision des couleurs et le contraste sur fond noir.
const COULEURS_ROI = [
    colorant"#3987e5", colorant"#d95926", colorant"#199e70", colorant"#c98500",
    colorant"#d55181", colorant"#008300", colorant"#9085e9", colorant"#e66767",
]

# --- Graphiques synchronises -------------------------------------------------
# Panneau de droite, meme taille que l'image du microscope : position x des
# galvos, puissance des deux lasers, puis les voies numeriques du generateur :
# experience active, impulsion de debut de sequence, impulsion de debut de ROI,
# laser 1 actif, laser 2 actif, aucun laser actif, et le numero de ROI code sur
# quatre bits. Le tout empile verticalement.
const AFFICHER_PLOTS   = true    # false = image du microscope seule
const FENETRE_PLOTS_MS = 0.0     # 0 = toute la sequence, tracee au fur et a
                                 # mesure derriere une tete de lecture ;
                                 # sinon largeur de la fenetre glissante [ms]
const COULEUR_PLOT_X = colorant"#B8C2D1"   # trace de la position x
const COULEUR_PLOT_Y = colorant"#B8C2D1"   # trace de la position y
const COULEUR_PLOT_NUM = colorant"#B8C2D1"  # voies numeriques generiques
const COULEUR_PLOT_ROI = colorant"#4fc3f7"  # compteur binaire de ROI

# Largeur des impulsions de debut de sequence et de debut de ROI [ms].
const DUREE_IMPULSION_MS = 2.0

# =============================================================================
#  REGLAGES AVANCES
# =============================================================================

const LUMINOSITE       = 1.0               # gain avant compression tonale ;
                                           # monter pour saturer vers le blanc
const FOND             = colorant"#000000" # couleur de fond
const INTERPOLER_TRAJET = false            # true : deplacement continu entre
                                           # les points au lieu de sauts
const AFFICHER_TEMPS   = true              # compteur de temps en haut a gauche
const AFFICHER_ECHELLE = true              # barre d'echelle en bas a droite
const COULEUR_TEXTE    = RGB(0.72, 0.76, 0.82)
const TAILLE_TEXTE     = 4                 # police 5x7 de la surcouche sur
                                           # l'image du microscope
const TAILLE_TEXTE_PLOTS = 3               # police des noms de voies, plus
                                           # petite car les bandes sont serrees
const RAYON_SPOT_SIGMA = 4.0               # rayon de calcul du spot, en sigmas
const SUR_ECHANTILLONNAGE = 4.0            # sous-pas par colonne de trace, pour
                                           # ne pas rater une impulsion breve
const QUALITE_H264     = 18                # crf ffmpeg ; plus bas = meilleur
const COULEUR_AXES     = RGB(0.34, 0.37, 0.42)  # cadres et lignes de zero
const EPAISSEUR_TRACE  = 0                 # epaisseur des traces en px ;
                                           # 0 = choisie selon la resolution.
                                           # Les axes et la tete de lecture
                                           # suivent proportionnellement.

# =============================================================================
#  (fin des parametres)
# =============================================================================

# -----------------------------------------------------------------------------
#  Trajectoire
# -----------------------------------------------------------------------------

"""
Points d'une tornade : spirale d'Archimede simple, non compensee, de rayon nul
au centre jusqu'a `rayon` au bord, parcourant `tours` tours complets.
`sortante` donne le sens radial : du centre vers le bord, ou l'inverse. La
rotation, elle, reste anti-horaire dans les deux cas ; seul le rayon change de
sens, l'angle croit toujours.
"""
function points_tornade(centre, n_points::Int, rayon::Real, tours::Real,
                        sortante::Bool)
    pts = Vector{NTuple{2,Float64}}(undef, n_points)
    for k in 1:n_points
        u = n_points == 1 ? 0.0 : (k - 1) / (n_points - 1)
        rho = rayon * (sortante ? u : 1.0 - u)
        theta = 2pi * tours * u              # toujours croissant : anti-horaire
        pts[k] = (centre[1] + rho * cos(theta), centre[2] + rho * sin(theta))
    end
    return pts
end

"""
Position des galvos au temps reel `t` [ms], et laser actif.
Retourne `(x, y, laser, roi)` avec `laser` valant 0 (eteint), 1 ou 2. Chaque laser a
sa propre trajectoire, les deux spirales pouvant tourner en sens opposes.
Pendant le delai les galvos sont deja sur le premier point de la prochaine
tornade, celui du laser qui ouvre la marche, les deux lasers eteints.
"""
function etat_faisceau(t::Float64, traj1, traj2, t_unite, dwell)
    n_tor = length(traj1)
    ouvre = INVERSER_ORDRE_LASERS ? 2 : 1          # laser de la premiere passe
    ferme = INVERSER_ORDRE_LASERS ? 1 : 2

    if t < 0.0
        p = (ouvre == 1 ? traj1 : traj2)[1][1]
        return (p[1], p[2], 0, 1)
    end

    creneau = floor(Int, t / t_unite)
    t_local = t - creneau * t_unite

    roi = mod(creneau, n_tor) + 1

    if t_local >= 2 * DUREE_TORNADE_MS              # delai : lasers eteints
        suiv = (ouvre == 1 ? traj1 : traj2)[mod(creneau + 1, n_tor)+1]
        return (suiv[1][1], suiv[1][2], 0, roi)
    end

    # deux passes enchainees sans pause, chacune avec son propre sens
    premiere = t_local < DUREE_TORNADE_MS
    laser = premiere ? ouvre : ferme
    t_passe = premiere ? t_local : t_local - DUREE_TORNADE_MS
    pts = (laser == 1 ? traj1 : traj2)[roi]

    if INTERPOLER_TRAJET
        s = clamp(t_passe / DUREE_TORNADE_MS, 0.0, 1.0) * (N_POINTS - 1) + 1
        i = clamp(floor(Int, s), 1, N_POINTS - 1)
        f = s - i
        x = pts[i][1] + f * (pts[i+1][1] - pts[i][1])
        y = pts[i][2] + f * (pts[i+1][2] - pts[i][2])
        return (x, y, laser, roi)
    else
        i = clamp(floor(Int, t_passe / dwell) + 1, 1, N_POINTS)
        return (pts[i][1], pts[i][2], laser, roi)
    end
end

# -----------------------------------------------------------------------------
#  Rendu
# -----------------------------------------------------------------------------

@inline function depose_spot!(H, xp::Float64, yp::Float64, sigma_px::Float64,
                              poids::Float64, rayon_px::Float64)
    W, Hh = size(H)
    i0 = max(1, ceil(Int, xp - rayon_px));  i1 = min(W,  floor(Int, xp + rayon_px))
    j0 = max(1, ceil(Int, yp - rayon_px));  j1 = min(Hh, floor(Int, yp + rayon_px))
    (i0 > i1 || j0 > j1) && return
    inv2s2 = 1.0 / (2 * sigma_px^2)
    r2max = rayon_px^2
    @inbounds for j in j0:j1
        dy = j - yp
        dy2 = dy * dy
        for i in i0:i1
            dx = i - xp
            r2 = dx * dx + dy2
            r2 > r2max && continue
            H[i, j] += poids * exp(-r2 * inv2s2)
        end
    end
    return
end

"""
Compression tonale puis coloration. Chaque laser garde sa teinte propre et la
dose la plus forte atteint exactement la couleur demandee, sans deborder vers
le blanc. Les deux contributions s'additionnent la ou les remanences se
recouvrent, comme deux sources lumineuses.
"""
function colorise!(img, H1, H2, c1, c2, fond)
    ar, ag, ab = Float64(red(c1)), Float64(green(c1)), Float64(blue(c1))
    br, bg, bb = Float64(red(c2)), Float64(green(c2)), Float64(blue(c2))
    fr, fg, fb = Float64(red(fond)), Float64(green(fond)), Float64(blue(fond))
    W, Hh = size(H1)
    @inbounds for j in 1:Hh, i in 1:W
        v1 = 1.0 - exp(-LUMINOSITE * H1[i, j])
        v2 = 1.0 - exp(-LUMINOSITE * H2[i, j])
        r = fr + (ar - fr) * v1 + (br - fr) * v2
        g = fg + (ag - fg) * v1 + (bg - fg) * v2
        b = fb + (ab - fb) * v1 + (bb - fb) * v2
        img[1, i, j] = round(UInt8, 255 * clamp(r, 0.0, 1.0))
        img[2, i, j] = round(UInt8, 255 * clamp(g, 0.0, 1.0))
        img[3, i, j] = round(UInt8, 255 * clamp(b, 0.0, 1.0))
    end
    return
end

# -----------------------------------------------------------------------------
#  Surcouche : police 5x7 minimale, barre d'echelle
# -----------------------------------------------------------------------------

const FONTE = Dict{Char,NTuple{7,UInt8}}(
    '0' => (0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110),
    '1' => (0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110),
    '2' => (0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111),
    '3' => (0b11111, 0b00010, 0b00100, 0b00010, 0b00001, 0b10001, 0b01110),
    '4' => (0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010),
    '5' => (0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110),
    '6' => (0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110),
    '7' => (0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000),
    '8' => (0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110),
    '9' => (0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100),
    '.' => (0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b01100, 0b01100),
    ':' => (0b00000, 0b01100, 0b01100, 0b00000, 0b01100, 0b01100, 0b00000),
    '=' => (0b00000, 0b00000, 0b11111, 0b00000, 0b11111, 0b00000, 0b00000),
    '-' => (0b00000, 0b00000, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000),
    '/' => (0b00001, 0b00010, 0b00010, 0b00100, 0b01000, 0b01000, 0b10000),
    ' ' => (0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000),
    'm' => (0b00000, 0b00000, 0b11010, 0b10101, 0b10101, 0b10101, 0b10101),
    's' => (0b00000, 0b00000, 0b01111, 0b10000, 0b01110, 0b00001, 0b11110),
    't' => (0b01000, 0b01000, 0b11110, 0b01000, 0b01000, 0b01001, 0b00110),
    'V' => (0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100),
    'a' => (0b00000, 0b00000, 0b01110, 0b00001, 0b01111, 0b10001, 0b01111),
    'e' => (0b00000, 0b00000, 0b01110, 0b10001, 0b11111, 0b10000, 0b01110),
    'i' => (0b00100, 0b00000, 0b01100, 0b00100, 0b00100, 0b00100, 0b01110),
    'l' => (0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110),
    'n' => (0b00000, 0b00000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001),
    'o' => (0b00000, 0b00000, 0b01110, 0b10001, 0b10001, 0b10001, 0b01110),
    'p' => (0b00000, 0b00000, 0b11110, 0b10001, 0b11110, 0b10000, 0b10000),
    'r' => (0b00000, 0b00000, 0b10110, 0b11001, 0b10000, 0b10000, 0b10000),
    'u' => (0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b10011, 0b01101),
    'v' => (0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100),
    'x' => (0b00000, 0b00000, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001),
    'y' => (0b00000, 0b00000, 0b10001, 0b10001, 0b01111, 0b00001, 0b01110),
    'I' => (0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110),
    'M' => (0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001),
    'O' => (0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110),
    'R' => (0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001),
    'b' => (0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b10001, 0b11110),
    'c' => (0b00000, 0b00000, 0b01111, 0b10000, 0b10000, 0b10000, 0b01111),
    'd' => (0b00001, 0b00001, 0b01111, 0b10001, 0b10001, 0b10001, 0b01111),
    'f' => (0b00110, 0b01001, 0b01000, 0b11100, 0b01000, 0b01000, 0b01000),
    'g' => (0b00000, 0b00000, 0b01111, 0b10001, 0b01111, 0b00001, 0b01110),
    'j' => (0b00010, 0b00000, 0b00110, 0b00010, 0b00010, 0b10010, 0b01100),
    'h' => (0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b10001, 0b10001),
    'k' => (0b10000, 0b10001, 0b10010, 0b11100, 0b10010, 0b10001, 0b10001),
    'q' => (0b00000, 0b00000, 0b01111, 0b10001, 0b01111, 0b00001, 0b00001),
    'w' => (0b00000, 0b00000, 0b10001, 0b10001, 0b10101, 0b10101, 0b01010),
    'z' => (0b00000, 0b00000, 0b11111, 0b00010, 0b00100, 0b01000, 0b11111),
    '(' => (0b00010, 0b00100, 0b01000, 0b01000, 0b01000, 0b00100, 0b00010),
    ')' => (0b01000, 0b00100, 0b00010, 0b00010, 0b00010, 0b00100, 0b01000),
)

@inline function peint!(img, i, j, c)
    W, Hh = size(img, 2), size(img, 3)
    (1 <= i <= W && 1 <= j <= Hh) || return
    @inbounds img[1, i, j] = c[1]
    @inbounds img[2, i, j] = c[2]
    @inbounds img[3, i, j] = c[3]
    return
end

function ecris_texte!(img, texte::AbstractString, x0::Int, y0::Int, ech::Int, c)
    x = x0
    for ch in texte
        glyphe = get(FONTE, ch, FONTE[' '])
        for (ligne, bits) in enumerate(glyphe), col in 1:5
            (bits >> (5 - col)) & 0x01 == 0x01 || continue
            for dy in 0:ech-1, dx in 0:ech-1
                peint!(img, x + (col - 1) * ech + dx, y0 + (ligne - 1) * ech + dy, c)
            end
        end
        x += 6 * ech
    end
    return x
end

largeur_texte(texte, ech) = 6 * ech * length(texte)

"Longueur ronde (1, 2 ou 5 x 10^k) valant environ 20 % de `span`."
function longueur_echelle(span::Float64)
    brut = 0.2 * span
    p = 10.0^floor(log10(brut))
    for m in (1.0, 2.0, 5.0)
        m * p >= brut && return m * p
    end
    return 10.0 * p
end

function barre_echelle!(img, W::Int, Hh::Int, ech_px::Float64, c)
    L_mV = longueur_echelle(XLIM[2] - XLIM[1])
    L_px = round(Int, L_mV * ech_px)
    marge = max(12, round(Int, 0.025 * W))
    ep = max(2, TAILLE_TEXTE ÷ 2)
    x1 = W - marge
    x0 = x1 - L_px
    y = Hh - marge - 7 * TAILLE_TEXTE - 5 * ep

    for i in x0:x1, dy in 0:ep-1
        peint!(img, i, y + dy, c)
    end
    for j in (y - 3 * ep):(y + 3 * ep), dx in 0:ep-1
        peint!(img, x0 + dx, j, c)
        peint!(img, x1 - dx, j, c)
    end

    label = string(round(Int, L_mV), " mV")
    ecris_texte!(img, label, x0 + (L_px - largeur_texte(label, TAILLE_TEXTE)) ÷ 2,
                 y + 3 * ep + 2 * ep, TAILLE_TEXTE, c)
    return
end

# -----------------------------------------------------------------------------
#  Graphiques synchronises (panneau de droite)
#
#  Le panneau reproduit les sorties d'un generateur de fonctions :
#  quatre voies analogiques puis dix voies numeriques, empilees. Les voies
#  analogiques recoivent le double de hauteur, les numeriques etant binaires.
# -----------------------------------------------------------------------------

# Poids de hauteur de chaque bande : 4 analogiques puis 10 numeriques
const POIDS_BANDES = (2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
const N_SIGNAUX = length(POIDS_BANDES)

octets(c) = (round(UInt8, 255 * red(c)),
             round(UInt8, 255 * green(c)),
             round(UInt8, 255 * blue(c)))

function rect!(img, x0::Int, x1::Int, y0::Int, y1::Int, c)
    for j in y0:y1, i in x0:x1
        peint!(img, i, j, c)
    end
    return
end

"Segment vertical epais, `x` au centre."
@inline function trait_v!(img, x::Int, y0::Int, y1::Int, c, ep::Int)
    d = ep ÷ 2
    rect!(img, x - d, x - d + ep - 1, min(y0, y1), max(y0, y1), c)
    return
end

"Segment horizontal epais, `y` au centre."
@inline function trait_h!(img, x0::Int, x1::Int, y::Int, c, ep::Int)
    d = ep ÷ 2
    rect!(img, min(x0, x1), max(x0, x1), y - d, y - d + ep - 1, c)
    return
end

"""
Etat complet du generateur au temps reel `t` [ms], dans l'ordre d'empilement :

  1-2   galvo x, galvo y                    (analogique, mV)
  3     puissance laser 1                    (analogique, 0 ou 1)
  4     puissance laser 2                    (analogique, commande du PI)
  5     experience demarree, puis maintenu  (numerique)
  6     impulsion au debut d'une sequence   (numerique)
  7     impulsion au debut d'une ROI        (numerique)
  8-10  laser 1 actif, laser 2 actif, aucun laser actif
  11-14 numero de ROI en binaire, bit 0 en tete, valide seulement quand un
        laser est actif. La ROI 1 vaut 1, donc 0001.
"""
function signaux(t::Float64, traj1, traj2, t_unite, dwell, n_tor::Int,
                 duree_boucle::Float64, xp)
    x, y, laser, roi = etat_faisceau(t, traj1, traj2, t_unite, dwell)
    actif = laser != 0

    creneau = t < 0.0 ? 0 : floor(Int, t / t_unite)
    t_roi = t - creneau * t_unite
    u2 = laser == 2 ? commande_a(xp, roi, t) : 0.0
    t_seq = t - floor(t / duree_boucle) * duree_boucle

    demarre = t >= 0.0
    imp_seq = (demarre && t_seq < DUREE_IMPULSION_MS) ? 1.0 : 0.0
    imp_roi = (demarre && t_roi < DUREE_IMPULSION_MS) ? 1.0 : 0.0
    bit(b) = (actif && ((roi >> (b - 1)) & 1) == 1) ? 1.0 : 0.0

    return (x, y,
            laser == 1 ? 1.0 : 0.0,
            u2,                                    # amplitude issue du PI
            demarre ? 1.0 : 0.0,
            imp_seq,
            imp_roi,
            laser == 1 ? 1.0 : 0.0,
            laser == 2 ? 1.0 : 0.0,
            laser == 0 ? 1.0 : 0.0,
            bit(1), bit(2), bit(3), bit(4))
end

"""
Nom, echelle et style de chaque voie, dans l'ordre d'empilement. Le detail
n'est rempli que pour les voies analogiques, ou l'echelle n'est pas evidente.
"""
function descripteurs_plots()
    ex = @sprintf("%.0f a %.0f mV", XLIM[1], XLIM[2])
    ey = @sprintf("%.0f a %.0f mV", YLIM[1], YLIM[2])
    num(nom, coul) = (nom = nom, detail = "", vmin = -0.12, vmax = 1.12,
                      coul = coul, remplir = true, zero = false)
    return [
        (nom = "galvo x", detail = ex, vmin = XLIM[1], vmax = XLIM[2],
         coul = COULEUR_PLOT_X, remplir = false, zero = true),
        (nom = "galvo y", detail = ey, vmin = YLIM[1], vmax = YLIM[2],
         coul = COULEUR_PLOT_Y, remplir = false, zero = true),
        (nom = "puissance 1", detail = "0 a 1", vmin = -0.12, vmax = 1.12,
         coul = COULEUR_LASER_1, remplir = true, zero = false),
        (nom = "puissance 2", detail = "0 a 1", vmin = -0.12, vmax = 1.12,
         coul = COULEUR_LASER_2, remplir = true, zero = false),
        num("exp. active", COULEUR_PLOT_NUM),
        num("debut sequence", COULEUR_PLOT_NUM),
        num("debut ROI", COULEUR_PLOT_NUM),
        num("laser 1 actif", COULEUR_LASER_1),
        num("laser 2 actif", COULEUR_LASER_2),
        num("aucun laser", COULEUR_PLOT_NUM),
        num("ROI bit 0", COULEUR_PLOT_ROI),
        num("ROI bit 1", COULEUR_PLOT_ROI),
        num("ROI bit 2", COULEUR_PLOT_ROI),
        num("ROI bit 3", COULEUR_PLOT_ROI),
    ]
end

"""
Decoupe le panneau de droite en bandes ponderees par POIDS_BANDES et retourne,
pour chacune, le rectangle de trace `(px0, px1, py0, py1)` et le haut de la
bande. Une reserve est gardee en bas pour l'axe du temps.
"""
function zones_plots(W_img::Int, H_img::Int)
    T = TAILLE_TEXTE_PLOTS
    # la marge gauche doit loger le plus long libelle, quelle que soit ECHELLE_PX
    plus_long = maximum(max(largeur_texte(v.nom, T), largeur_texte(v.detail, T))
                        for v in descripteurs_plots())
    px0 = W_img + max(round(Int, 0.24 * W_img), plus_long + 8 * T)
    px1 = 2 * W_img - round(Int, 0.03 * W_img)
    utile = H_img - 10 * T                      # reserve pour l'axe du temps
    total = sum(POIDS_BANDES)
    zones = Vector{NTuple{5,Int}}(undef, N_SIGNAUX)
    cumul = 0
    for p in 1:N_SIGNAUX
        haut = round(Int, cumul / total * utile) + 1
        cumul += POIDS_BANDES[p]
        bas = round(Int, cumul / total * utile)
        zones[p] = (px0, px1, haut + 3 * T, bas - 2 * T, haut)
    end
    return zones
end

"""
Pour chaque colonne de pixel et chaque voie, l'intervalle `(min, max)` pris par
le signal pendant l'intervalle de temps que cette colonne recouvre.
Echantillonner plusieurs fois par colonne evite le crenelage quand la spirale
oscille vite, et garantit qu'une impulsion breve n'est pas ratee.

"""
function colonnes_traces(t_a::Float64, t_b::Float64, n_col::Int,
                         traj1, traj2, t_unite, dwell, n_tor::Int,
                         duree_boucle::Float64, xp, n_ech::Int)
    series = [Vector{NTuple{2,Float64}}(undef, n_col) for _ in 1:N_SIGNAUX]
    dt = (t_b - t_a) / n_col
    lo = fill(Inf, N_SIGNAUX)
    hi = fill(-Inf, N_SIGNAUX)
    for c in 1:n_col
        fill!(lo, Inf)
        fill!(hi, -Inf)
        ta = t_a + (c - 1) * dt
        for e in 1:n_ech
            t = ta + (e - 0.5) * dt / n_ech
            v = signaux(t, traj1, traj2, t_unite, dwell, n_tor, duree_boucle, xp)
            @inbounds for i in 1:N_SIGNAUX
                lo[i] = min(lo[i], v[i])
                hi[i] = max(hi[i], v[i])
            end
        end
        @inbounds for i in 1:N_SIGNAUX
            series[i][c] = (lo[i], hi[i])
        end
    end
    return series
end

"""
Une voie : cadre, ligne de zero, trace de toute la sequence en cours, et nom de
la voie dans la marge de gauche. Le trace est peint d'un seul tenant, sans tete
de lecture ni fondu.
"""
function dessine_plot!(img, zone, donnees, d, ep::Int, ep_axe::Int, axes_u8)
    px0, px1, py0, py1, haut = zone
    T = TAILLE_TEXTE_PLOTS
    coul = d.coul
    c_u8 = octets(coul)
    sombre = let cr = Float64(red(coul)), cg = Float64(green(coul)),
                 cb = Float64(blue(coul)),
                 fr = Float64(red(FOND)), fg = Float64(green(FOND)),
                 fb = Float64(blue(FOND))
        (round(UInt8, 255 * (fr + (cr - fr) * 0.30)),
         round(UInt8, 255 * (fg + (cg - fg) * 0.30)),
         round(UInt8, 255 * (fb + (cb - fb) * 0.30)))
    end

    rangee(v) = clamp(round(Int, py1 - (v - d.vmin) / (d.vmax - d.vmin) *
                                  (py1 - py0)), py0, py1)

    # cadre
    trait_h!(img, px0, px1, py0, axes_u8, ep_axe)
    trait_h!(img, px0, px1, py1, axes_u8, ep_axe)
    trait_v!(img, px0, py0, py1, axes_u8, ep_axe)

    if d.zero
        rz = rangee(0.0)
        for i in px0:(4 * T):px1                      # pointilles
            trait_h!(img, i, min(i + 2 * T, px1), rz, axes_u8, ep_axe)
        end
    end

    base = rangee(0.0)
    dd = ep ÷ 2
    n_tot = length(donnees)

    # rangees brutes de chaque colonne, puis pontage avec la precedente pour
    # que les fronts raides soient des segments continus
    blo = Vector{Int}(undef, n_tot)
    bhi = Vector{Int}(undef, n_tot)
    for c in 1:n_tot
        lo_v, hi_v = donnees[c]
        blo[c] = rangee(hi_v)               # valeur haute -> rangee basse
        bhi[c] = rangee(lo_v)
    end
    rlo, rhi = copy(blo), copy(bhi)
    for c in 2:n_tot
        rlo[c] = min(blo[c], bhi[c-1])
        rhi[c] = max(bhi[c], blo[c-1])
    end

    # deux passes : les remplissages d'abord, sinon le pinceau large d'une
    # colonne pleine efface le front qui la precede
    if d.remplir
        for c in 1:n_tot
            x = px0 + c - 1
            rect!(img, x, x, min(base, rlo[c]), max(base, rlo[c]), sombre)
        end
    end
    # pinceau carre : meme epaisseur sur les paliers et sur les fronts
    for c in 1:n_tot
        x = px0 + c - 1
        rect!(img, x - dd, x - dd + ep - 1,
              clamp(rlo[c] - dd, py0, py1),
              clamp(rhi[c] - dd + ep - 1, py0, py1), c_u8)
    end

    # nom de la voie, aligne a droite dans la marge de gauche
    cy = (py0 + py1) ÷ 2
    if isempty(d.detail)
        ecris_texte!(img, d.nom, px0 - largeur_texte(d.nom, T) - 3 * T,
                     cy - 3 * T, T, c_u8)
    else
        ecris_texte!(img, d.nom, px0 - largeur_texte(d.nom, T) - 3 * T,
                     cy - 8 * T, T, c_u8)
        ecris_texte!(img, d.detail, px0 - largeur_texte(d.detail, T) - 3 * T,
                     cy + 2 * T, T, axes_u8)
    end
    return
end

"Peint le panneau de droite : toutes les voies du generateur, empilees."
function dessine_plots!(img, W_img::Int, H_img::Int, zones, descr, traces,
                        ep::Int, t_a::Float64, t_b::Float64)
    fond_u8 = octets(FOND)
    axes_u8 = octets(COULEUR_AXES)
    ep_axe = max(2, ep ÷ 4)
    T = TAILLE_TEXTE_PLOTS
    rect!(img, W_img + 1, 2 * W_img, 1, H_img, fond_u8)
    trait_v!(img, W_img + 1, 1, H_img, axes_u8, ep_axe)

    for p in 1:N_SIGNAUX
        dessine_plot!(img, zones[p], traces[p], descr[p], ep, ep_axe, axes_u8)
    end

    # axe du temps sous la bande du bas
    px0, px1, _, py1, _ = zones[end]
    g = @sprintf("%.0f ms", t_a)
    dr = @sprintf("%.0f ms", t_b)
    ecris_texte!(img, g, px0, py1 + 3 * T, T, axes_u8)
    ecris_texte!(img, dr, px1 - largeur_texte(dr, T), py1 + 3 * T, T, axes_u8)
    return
end

# -----------------------------------------------------------------------------
#  Experience : acquisitions, correcteur PI et reponse du premier ordre
# -----------------------------------------------------------------------------

"""
Simule toute l'experience avant le rendu. Pour chaque ROI, retourne les instants
d'acquisition [ms], la concentration mesuree et la commande appliquee jusqu'a
l'acquisition suivante.

L'acquisition d'une ROI tombe a la fin de sa passe de laser 1 ; la commande
calculee est appliquee immediatement, donc des la passe de laser 2 qui suit.
Entre deux acquisitions la commande est bloquee, ce qui permet d'integrer le
premier ordre exactement au lieu de le discretiser.
"""
function simule_clampage(n_tor::Int, t_unite::Float64, t_seq::Float64, n_acq::Int)
    dt = t_seq / 1000.0                      # periode d'acquisition [s]
    rng = MersenneTwister(GRAINE_BRUIT)

    t_acq = [Vector{Float64}(undef, n_acq) for _ in 1:n_tor]
    cl    = [Vector{Float64}(undef, n_acq) for _ in 1:n_tor]
    cmd   = [Vector{Float64}(undef, n_acq) for _ in 1:n_tor]

    for i in 1:n_tor
        base = CL_BASAL[i]
        tau = TAU_ROI_S[i]
        etat = base                          # concentration vraie
        integ = 0.0
        decalage = (i - 1) * t_unite + DUREE_TORNADE_MS

        for k in 1:n_acq
            t_acq[i][k] = (k - 1) * t_seq + decalage
            mesure = etat + (BRUIT_MESURE > 0 ? BRUIT_MESURE * randn(rng) : 0.0)
            cl[i][k] = mesure

            if CLAMP_DEBUT_ACQ <= k < CLAMP_FIN_ACQ
                e = CL_CIBLE - mesure
                brut = PI_KP * e + PI_KI * integ
                u = clamp(brut, 0.0, 1.0)
                # Integration conditionnelle : on ne gele l'integrale que si la
                # commande sature ET que l'erreur pousse encore dans ce sens.
                # La geler des que la commande brute sort de [0, 1] la bloquerait
                # a zero tant que le terme proportionnel est grand, et le
                # correcteur ne rattraperait jamais l'ecart statique.
                pousse = (brut > 1.0 && e > 0.0) || (brut < 0.0 && e < 0.0)
                pousse || (integ += e * dt)
            else
                u = 0.0
                integ = 0.0
            end
            cmd[i][k] = u

            # reponse exacte du premier ordre, commande bloquee sur le pas
            cible = base + GAIN_COMMANDE * u
            etat = cible + (etat - cible) * exp(-dt / tau)
        end
    end
    return (t_acq = t_acq, cl = cl, cmd = cmd,
            t_seq = t_seq, t_unite = t_unite, n_acq = n_acq)
end

"Commande en vigueur pour la ROI `i` au temps `t` [ms], bloqueur d'ordre zero."
@inline function commande_a(xp, i::Int, t::Float64)
    decalage = (i - 1) * xp.t_unite + DUREE_TORNADE_MS
    k = floor(Int, (t - decalage) / xp.t_seq) + 1
    (k < 1 || k > xp.n_acq) && return 0.0
    return @inbounds xp.cmd[i][k]
end

# -----------------------------------------------------------------------------
#  Courbes du bas : concentration et commande, toutes les ROI superposees
# -----------------------------------------------------------------------------

"Decoupe la moitie basse en deux bandes et retourne leurs rectangles de trace."
function zones_courbes(W_tot::Int, H_img::Int, H_tot::Int)
    T = TAILLE_TEXTE_PLOTS
    # la marge gauche doit loger le plus long libelle, quelle que soit ECHELLE_PX
    plus_long = maximum(largeur_texte(m, T)
                        for m in ("concentration", "commande", "0 a 1"))
    px0 = max(round(Int, 0.13 * W_tot), plus_long + 10 * T)
    px1 = W_tot - round(Int, 0.03 * W_tot)
    utile = (H_tot - H_img) - 10 * T         # reserve pour l'axe du temps
    zones = Vector{NTuple{5,Int}}(undef, 2)
    for p in 1:2
        haut = H_img + round(Int, (p - 1) / 2 * utile) + 1
        bas  = H_img + round(Int, p / 2 * utile)
        zones[p] = (px0, px1, haut + 9 * T, bas - 9 * T, haut)
    end
    return zones
end

"""
Une courbe : cadre, une polyligne par ROI revelee jusqu'a `t_cur`, tete de
lecture, titre et legende dans la marge de gauche. Les marqueurs ne sont poses
que si les acquisitions sont assez espacees pour ne pas se confondre.
"""
function dessine_courbe!(img, zone, xp, valeurs, vmin::Float64, vmax::Float64,
                         titre::String, unite::String, etiq_max::String,
                         etiq_min::String, t_fin::Float64,
                         t_cur::Float64, ep::Int, ep_axe::Int, axes_u8,
                         legende::Bool)
    px0, px1, py0, py1, haut = zone
    T = TAILLE_TEXTE_PLOTS
    n_tor = length(valeurs)

    colonne(t) = px0 + round(Int, clamp(t / t_fin, 0.0, 1.0) * (px1 - px0))
    rangee(v) = clamp(round(Int, py1 - (v - vmin) / (vmax - vmin) * (py1 - py0)),
                      py0, py1)

    trait_h!(img, px0, px1, py0, axes_u8, ep_axe)
    trait_h!(img, px0, px1, py1, axes_u8, ep_axe)
    trait_v!(img, px0, py0, py1, axes_u8, ep_axe)

    # marqueurs seulement si les points ne se touchent pas
    pas_px = (px1 - px0) / max(1, xp.n_acq - 1)
    marqueur = pas_px >= 3.5 * ep ? max(ep + 2, 8) : 0

    for i in 1:n_tor
        c_u8 = octets(COULEURS_ROI[i])
        xs = xp.t_acq[i]
        vs = valeurs[i]
        xp_prec = 0
        yp_prec = 0
        premier = true
        for k in 1:xp.n_acq
            xs[k] > t_cur && break
            xk = colonne(xs[k])
            yk = rangee(vs[k])
            if premier
                premier = false
            else
                segment!(img, xp_prec, yp_prec, xk, yk, c_u8, ep)
            end
            marqueur > 0 && rect!(img, xk - marqueur ÷ 2, xk + marqueur ÷ 2,
                                  yk - marqueur ÷ 2, yk + marqueur ÷ 2, c_u8)
            xp_prec, yp_prec = xk, yk
        end
    end

    # tete de lecture
    xc = colonne(t_cur)
    px0 < xc < px1 && trait_v!(img, xc, py0, py1, axes_u8, max(2, ep ÷ 2))

    # marge de gauche, de haut en bas : borne haute, titre, unite, legende,
    # borne basse. Les libelles sont alignes a droite sur le bord du cadre.
    droite(txt) = px0 - largeur_texte(txt, T) - 4 * T
    haut_txt, bas_txt = etiq_max, etiq_min
    ecris_texte!(img, haut_txt, droite(haut_txt), py0 - 8 * T, T, axes_u8)
    ecris_texte!(img, titre, droite(titre), py0 + 3 * T, T, octets(COULEUR_TEXTE))
    ecris_texte!(img, unite, droite(unite), py0 + 12 * T, T, axes_u8)
    ecris_texte!(img, bas_txt, droite(bas_txt), py1 - 8 * T, T, axes_u8)
    if legende
        for i in 1:n_tor
            y = py0 + (24 + 10 * (i - 1)) * T
            nom = "ROI " * string(i)
            xg = droite(nom)
            rect!(img, xg - 7 * T, xg - 3 * T, y, y + 5 * T,
                  octets(COULEURS_ROI[i]))
            ecris_texte!(img, nom, xg, y - T, T, octets(COULEUR_TEXTE))
        end
    end
    return
end

"Segment epais entre deux points, trace colonne par colonne."
function segment!(img, x0::Int, y0::Int, x1::Int, y1::Int, c, ep::Int)
    d = ep ÷ 2
    if x1 == x0
        rect!(img, x0 - d, x0 - d + ep - 1, min(y0, y1) - d,
              max(y0, y1) - d + ep - 1, c)
        return
    end
    pente = (y1 - y0) / (x1 - x0)
    ya = y0
    for x in x0:x1
        yb = y0 + round(Int, pente * (x - x0))
        rect!(img, x - d, x - d + ep - 1, min(ya, yb) - d,
              max(ya, yb) - d + ep - 1, c)
        ya = yb
    end
    return
end

"Peint la moitie basse : concentration puis commande, toutes les ROI ensemble."
function dessine_courbes!(img, W_tot::Int, H_img::Int, H_tot::Int, zones, xp,
                          t_fin::Float64, t_cur::Float64, ep::Int)
    fond_u8 = octets(FOND)
    axes_u8 = octets(COULEUR_AXES)
    ep_axe = max(2, ep ÷ 4)
    T = TAILLE_TEXTE_PLOTS
    rect!(img, 1, W_tot, H_img + 1, H_tot, fond_u8)
    trait_h!(img, 1, W_tot, H_img + 1, axes_u8, ep_axe)

    dessine_courbe!(img, zones[1], xp, xp.cl, CL_PLOT[1], CL_PLOT[2],
                    "concentration", "mM",
                    @sprintf("%.0f", CL_PLOT[2]), @sprintf("%.0f", CL_PLOT[1]),
                    t_fin, t_cur, ep, ep_axe, axes_u8, true)
    dessine_courbe!(img, zones[2], xp, xp.cmd, -0.05, 1.05,
                    "commande", "0 a 1", "1", "0",
                    t_fin, t_cur, ep, ep_axe, axes_u8, false)

    # consigne du clampage, en pointilles sur le graphe de concentration
    px0, px1, py0, py1, _ = zones[1]
    rz = clamp(round(Int, py1 - (CL_CIBLE - CL_PLOT[1]) /
                           (CL_PLOT[2] - CL_PLOT[1]) * (py1 - py0)), py0, py1)
    for i in px0:(6 * T):px1
        trait_h!(img, i, min(i + 3 * T, px1), rz, axes_u8, ep_axe)
    end
    ecris_texte!(img, "consigne", px1 - largeur_texte("consigne", T),
                 rz - 9 * T, T, axes_u8)

    # axe du temps sous la bande du bas
    px0b, px1b, _, py1b, _ = zones[2]
    dr = @sprintf("%.1f s", t_fin / 1000)
    ecris_texte!(img, "0 s", px0b, py1b + 3 * T, T, axes_u8)
    ecris_texte!(img, dr, px1b - largeur_texte(dr, T), py1b + 3 * T, T, axes_u8)
    return
end

# -----------------------------------------------------------------------------
#  Simulation
# -----------------------------------------------------------------------------

function simule()
    isempty(POSITIONS) && error("POSITIONS est vide : aucune tornade a simuler.")
    n_pos = length(POSITIONS)
    length(CL_BASAL) == n_pos || error(
        "CL_BASAL a $(length(CL_BASAL)) valeurs pour $n_pos positions.")
    length(TAU_ROI_S) == n_pos || error(
        "TAU_ROI_S a $(length(TAU_ROI_S)) valeurs pour $n_pos positions.")
    n_pos > length(COULEURS_ROI) && error(
        "La palette des courbes compte $(length(COULEURS_ROI)) couleurs, " *
        "jamais recyclees : $n_pos positions demandees.")
    length(POSITIONS) > 15 && error(
        "Le compteur de ROI est code sur 4 bits : 15 positions au maximum, " *
        "$(length(POSITIONS)) demandees.")
    N_POINTS < 1 && error("N_POINTS doit valoir au moins 1.")

    n_tor = length(POSITIONS)
    traj1 = [points_tornade(p, N_POINTS, RAYON, N_TOURS, SPIRALE_SORTANTE_1)
             for p in POSITIONS]
    traj2 = [points_tornade(p, N_POINTS, RAYON, N_TOURS, SPIRALE_SORTANTE_2)
             for p in POSITIONS]

    t_unite = 2 * DUREE_TORNADE_MS + DELAI_MS      # deux passes + le delai [ms]
    dwell = DUREE_TORNADE_MS / N_POINTS            # temps par point [ms]
    duree_reelle_ms = N_REPETITIONS * n_tor * t_unite
    duree_video_s = duree_reelle_ms / 1000 * RALENTI
    n_frames = max(1, round(Int, duree_video_s * FPS))
    dt_reel = duree_reelle_ms / n_frames           # temps reel par image [ms]

    W = round(Int, (XLIM[2] - XLIM[1]) * ECHELLE_PX) + 1
    Hpx = round(Int, (YLIM[2] - YLIM[1]) * ECHELLE_PX) + 1
    W_tot = AFFICHER_PLOTS ? 2W : W
    H_tot = 2Hpx                                   # moitie basse : les courbes

    sigma_px = max(0.6, SIGMA_FAISCEAU * ECHELLE_PX)
    rayon_px = RAYON_SPOT_SIGMA * sigma_px

    mkpath(DOSSIER_SORTIE)
    nom = isempty(NOM_FICHIER) ?
          "tornades_" * Dates.format(now(), raw"yyyy-mm-dd_HH\hMM\mSS") * ".mp4" :
          NOM_FICHIER
    chemin = joinpath(DOSSIER_SORTIE, nom)

    println("Tornades          : $n_tor position(s) x $N_REPETITIONS repetition(s)")
    println("Sequence          : laser 1, laser 2, delai, a chaque position")
    println("Voies generateur  : 4 analogiques + 10 numeriques")
    @printf("Acquisitions      : %d par ROI, clampage des la %d, relache a la %d\n",
            N_REPETITIONS, CLAMP_DEBUT_ACQ, CLAMP_FIN_ACQ)
    @printf("Temps reel        : %.0f ms (%.2f ms par point)\n", duree_reelle_ms, dwell)
    @printf("Video             : %d x %d px, %.0f fps, %.1f s (ralenti %gx)\n",
            W_tot, H_tot, FPS, n_frames / FPS, RALENTI)
    println("Images            : $n_frames")
    println("Sortie            : $chemin")
    if isodd(W_tot) || isodd(H_tot)
        println("Note              : dimension impaire, ffmpeg ajoute 1 px de bord noir.")
    end
    println()

    H1 = zeros(Float32, W, Hpx)                    # remanence du laser 1
    H2 = zeros(Float32, W, Hpx)                    # remanence du laser 2
    img = Array{UInt8}(undef, 3, W_tot, H_tot)
    txt = (round(UInt8, 255 * red(COULEUR_TEXTE)),
           round(UInt8, 255 * green(COULEUR_TEXTE)),
           round(UInt8, 255 * blue(COULEUR_TEXTE)))
    en_ms = duree_reelle_ms < 1000.0

    x_px(x) = (x - XLIM[1]) * ECHELLE_PX + 1.0
    y_px(y) = (YLIM[2] - y) * ECHELLE_PX + 1.0

    # graphiques : geometrie, echantillonnage et trace fixe si pas de fenetre
    zones = zones_plots(W, Hpx)
    descr = descripteurs_plots()
    n_col = zones[1][2] - zones[1][1] + 1
    ep_trace = EPAISSEUR_TRACE > 0 ? EPAISSEUR_TRACE : max(4, round(Int, W / 250))
    defilant = FENETRE_PLOTS_MS > 0
    duree_boucle_ms = n_tor * t_unite        # les plots montrent une boucle
    xp = simule_clampage(n_tor, t_unite, duree_boucle_ms, N_REPETITIONS)
    zones_bas = zones_courbes(W_tot, Hpx, H_tot)
    n_ech = clamp(ceil(Int, SUR_ECHANTILLONNAGE *
                       ((defilant ? FENETRE_PLOTS_MS : duree_boucle_ms) /
                        n_col) / dwell), 2, 64)


    taille = "$(W_tot)x$(H_tot)"
    vf = "pad=ceil(iw/2)*2:ceil(ih/2)*2:0:0:black"

    FFMPEG_jll.ffmpeg() do exe
        cmd = `$exe -hide_banner -loglevel error -y
               -f rawvideo -pix_fmt rgb24 -s $taille -r $FPS -i pipe:0
               -an -vf $vf -c:v libx264 -preset medium -crf $QUALITE_H264
               -pix_fmt yuv420p -movflags +faststart $chemin`

        open(cmd, "w") do flux
            for k in ProgressBar(1:n_frames)
                t1 = k * dt_reel                      # fin de l'image, temps reel

                # La tornade apparait et disparait d'un bloc : on peint tous
                # les points de la passe active a cet instant, sans remanence.
                # L'image obtenue est la dose integree d'une passe entiere.
                fill!(H1, 0.0f0)
                fill!(H2, 0.0f0)
                _, _, laser, roi = etat_faisceau(t1, traj1, traj2, t_unite, dwell)
                if laser != 0
                    poids = laser == 2 ? commande_a(xp, roi, t1) : 1.0
                    if poids > 0.0
                        pts = (laser == 1 ? traj1 : traj2)[roi]
                        acc = laser == 1 ? H1 : H2
                        for (bx, by) in pts
                            depose_spot!(acc, x_px(bx), y_px(by), sigma_px,
                                         poids, rayon_px)
                        end
                    end
                end

                colorise!(img, H1, H2, COULEUR_LASER_1, COULEUR_LASER_2, FOND)

                if AFFICHER_TEMPS
                    etiq = en_ms ? @sprintf("t = %6.1f ms", t1) :
                                   @sprintf("t = %6.3f s", t1 / 1000)
                    marge = max(12, round(Int, 0.025 * W))
                    ecris_texte!(img, etiq, marge, marge, TAILLE_TEXTE, txt)
                end
                AFFICHER_ECHELLE && barre_echelle!(img, W, Hpx, ECHELLE_PX, txt)

                if AFFICHER_PLOTS
                    if defilant
                        ta, tb = t1 - FENETRE_PLOTS_MS, t1
                        tr = colonnes_traces(ta, tb, n_col, traj1, traj2,
                                             t_unite, dwell, n_tor,
                                             duree_boucle_ms, xp, n_ech)
                        dessine_plots!(img, W, Hpx, zones, descr, tr,
                                       ep_trace, ta, tb)
                    else
                        # la sequence qui contient t1, peinte en entier
                        depart = floor(t1 / duree_boucle_ms) * duree_boucle_ms
                        tr = colonnes_traces(depart, depart + duree_boucle_ms,
                                             n_col, traj1, traj2, t_unite, dwell,
                                             n_tor, duree_boucle_ms, xp, n_ech)
                        dessine_plots!(img, W, Hpx, zones, descr, tr,
                                       ep_trace, 0.0, duree_boucle_ms)
                    end
                end

                dessine_courbes!(img, W_tot, Hpx, H_tot, zones_bas, xp,
                                 duree_reelle_ms, t1, ep_trace)

                write(flux, img)
            end
        end
    end

    println("\nVideo ecrite : $chemin")
    return chemin
end

simule()
