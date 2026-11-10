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
 gaussienne d'ecart-type SIGMA_FAISCEAU et la trainee est une remanence
 exponentielle de constante de temps TRAINEE_MS.

 Lancement : julia Simulateurs/Tornades.jl
===============================================================================
=#

using Colors
using Dates
using Printf
using ProgressBars
using FFMPEG_jll

# =============================================================================
#  PARAMETRES
# =============================================================================

# --- Geometrie d'une tornade -------------------------------------------------
const N_POINTS = 1000         # nombre de points par tornade
const N_TOURS  = 4.0          # nombre de tours complets de la spirale
                              # N_POINTS / N_TOURS donne les points par tour :
                              # sous environ 2*pi*RAYON/SIGMA_FAISCEAU points
                              # par tour, la spirale se lit en pointilles.
                              # INTERPOLER_TRAJET la rend continue.
const RAYON    = 100.0        # rayon exterieur de la tornade [mV]

# --- Sequence temporelle -----------------------------------------------------
const DUREE_TORNADE_MS = 200.0  # duree d'UNE passe [ms] ; chaque position en
                                # compte deux, enchainees sans pause
const DELAI_MS         = 100.0  # delai apres les deux passes, lasers eteints [ms]
const N_REPETITIONS    = 2      # nombre de repetitions de la boucle complete

# --- Positions des tornades --------------------------------------------------
# Le nombre de tornades est la longueur de cette liste. Elles sont visitees
# dans l'ordre ou elles sont ecrites.
const POSITIONS = [
    (-500.0,  500.0),
    (   0.0,  500.0),
    ( 500.0,  500.0),
    (-500.0,    0.0),
    (   0.0,    0.0),
    ( 500.0,    0.0),
    (-500.0, -500.0),
    (   0.0, -500.0),
    ( 500.0, -500.0),
]

# --- Faisceau ----------------------------------------------------------------
const SIGMA_FAISCEAU = 15.0                 # ecart-type du spot gaussien [mV]
const TRAINEE_MS     = 10.0                 # constante de temps de la trainee
                                            # [ms] ; 0.0 supprime la trainee
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
const FPS     = 30.0             # images par seconde de la video
const RALENTI = 0.1              # facteur de ralenti ; 1.0 = temps reel,
const DOSSIER_SORTIE = joinpath(@__DIR__, "sorties")
const NOM_FICHIER    = ""        # vide = nom horodate automatique

# --- Graphiques synchronises -------------------------------------------------
# Panneau de droite, meme taille que l'image du microscope : position x des
# galvos, position y des galvos, puissance du laser 1 et puissance du laser 2,
# empiles verticalement.
const AFFICHER_PLOTS   = true    # false = image du microscope seule
const FENETRE_PLOTS_MS = 0.0     # 0 = toute la sequence, tracee au fur et a
                                 # mesure derriere une tete de lecture ;
                                 # sinon largeur de la fenetre glissante [ms]
const COULEUR_PLOT_X = colorant"#B8C2D1"   # trace de la position x
const COULEUR_PLOT_Y = colorant"#B8C2D1"   # trace de la position y
const COULEUR_PLOT_P1 = COULEUR_LASER_1    # trace puissance laser 1
const COULEUR_PLOT_P2 = COULEUR_LASER_2    # trace puissance laser 2

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
const TAILLE_TEXTE     = 4                 # facteur d'agrandissement de la
                                           # police 5x7 de la surcouche
const RAYON_SPOT_SIGMA = 4.0               # rayon de calcul du spot, en sigmas
const SUR_ECHANTILLONNAGE = 2.0            # sous-pas de calcul par point visite
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
Retourne `(x, y, laser)` avec `laser` valant 0 (eteint), 1 ou 2. Chaque laser a
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
        return (p[1], p[2], 0)
    end

    creneau = floor(Int, t / t_unite)
    t_local = t - creneau * t_unite

    if t_local >= 2 * DUREE_TORNADE_MS              # delai : lasers eteints
        suiv = (ouvre == 1 ? traj1 : traj2)[mod(creneau + 1, n_tor)+1]
        return (suiv[1][1], suiv[1][2], 0)
    end

    # deux passes enchainees sans pause, chacune avec son propre sens
    premiere = t_local < DUREE_TORNADE_MS
    laser = premiere ? ouvre : ferme
    t_passe = premiere ? t_local : t_local - DUREE_TORNADE_MS
    pts = (laser == 1 ? traj1 : traj2)[mod(creneau, n_tor)+1]

    if INTERPOLER_TRAJET
        s = clamp(t_passe / DUREE_TORNADE_MS, 0.0, 1.0) * (N_POINTS - 1) + 1
        i = clamp(floor(Int, s), 1, N_POINTS - 1)
        f = s - i
        x = pts[i][1] + f * (pts[i+1][1] - pts[i][1])
        y = pts[i][2] + f * (pts[i+1][2] - pts[i][2])
        return (x, y, laser)
    else
        i = clamp(floor(Int, t_passe / dwell) + 1, 1, N_POINTS)
        return (pts[i][1], pts[i][2], laser)
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
    'b' => (0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b10001, 0b11110),
    'c' => (0b00000, 0b00000, 0b01111, 0b10000, 0b10000, 0b10000, 0b01111),
    'd' => (0b00001, 0b00001, 0b01111, 0b10001, 0b10001, 0b10001, 0b01111),
    'g' => (0b00000, 0b00000, 0b01111, 0b10001, 0b01111, 0b00001, 0b01110),
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

function barre_echelle!(img, W::Int, ech_px::Float64, c)
    Hh = size(img, 3)
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
# -----------------------------------------------------------------------------

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
Decoupe le panneau de droite en quatre bandes et retourne, pour chacune, le
rectangle de trace `(px0, px1, py0, py1)` et le haut de la bande.
"""
function zones_plots(W_img::Int, H_img::Int)
    T = TAILLE_TEXTE
    h_bande = H_img ÷ 4
    px0 = W_img + round(Int, 0.15 * W_img)
    px1 = 2 * W_img - round(Int, 0.04 * W_img)
    zones = Vector{NTuple{5,Int}}(undef, 4)
    for p in 1:4
        haut = (p - 1) * h_bande + 1
        zones[p] = (px0, px1, haut + 11 * T, haut + h_bande - 11 * T, haut)
    end
    return zones
end

"""
Pour chaque colonne de pixel du trace, l'intervalle `(min, max)` pris par le
signal pendant l'intervalle de temps que cette colonne recouvre. Echantillonner
plusieurs fois par colonne evite le crenelage quand la spirale oscille vite.
"""
function colonnes_traces(t_a::Float64, t_b::Float64, n_col::Int,
                         traj1, traj2, t_unite, dwell, n_ech::Int)
    gx = Vector{NTuple{2,Float64}}(undef, n_col)
    gy = Vector{NTuple{2,Float64}}(undef, n_col)
    g1 = Vector{NTuple{2,Float64}}(undef, n_col)
    g2 = Vector{NTuple{2,Float64}}(undef, n_col)
    dt = (t_b - t_a) / n_col
    for c in 1:n_col
        ta = t_a + (c - 1) * dt
        xlo = ylo = alo = blo = Inf
        xhi = yhi = ahi = bhi = -Inf
        for e in 1:n_ech
            t = ta + (e - 0.5) * dt / n_ech
            x, y, laser = etat_faisceau(t, traj1, traj2, t_unite, dwell)
            va = laser == 1 ? 1.0 : 0.0
            vb = laser == 2 ? 1.0 : 0.0
            xlo = min(xlo, x); xhi = max(xhi, x)
            ylo = min(ylo, y); yhi = max(yhi, y)
            alo = min(alo, va); ahi = max(ahi, va)
            blo = min(blo, vb); bhi = max(bhi, vb)
        end
        gx[c] = (xlo, xhi); gy[c] = (ylo, yhi)
        g1[c] = (alo, ahi); g2[c] = (blo, bhi)
    end
    return gx, gy, g1, g2
end

"""
Un graphique : cadre, ligne de zero, trace revele jusqu'a `c_rev` colonnes,
tete de lecture, titre et etiquettes d'axe.
"""
function dessine_plot!(img, zone, donnees, vmin::Float64, vmax::Float64,
                       coul, titre::String, etiq_haut::String, etiq_bas::String,
                       c_rev::Int, n_dessin::Int, ep::Int, ep_axe::Int,
                       remplir::Bool, axes_u8, ligne_zero::Bool)
    px0, px1, py0, py1, haut = zone
    T = TAILLE_TEXTE
    c_u8 = octets(coul)

    # teinte du trace : fondu vers le fond quand la colonne vieillit
    cr, cg, cb = Float64(red(coul)), Float64(green(coul)), Float64(blue(coul))
    fr, fg, fb = Float64(red(FOND)), Float64(green(FOND)), Float64(blue(FOND))
    teinte(a) = (round(UInt8, 255 * clamp(fr + (cr - fr) * a, 0.0, 1.0)),
                 round(UInt8, 255 * clamp(fg + (cg - fg) * a, 0.0, 1.0)),
                 round(UInt8, 255 * clamp(fb + (cb - fb) * a, 0.0, 1.0)))

    rangee(v) = clamp(round(Int, py1 - (v - vmin) / (vmax - vmin) * (py1 - py0)),
                      py0, py1)

    # cadre
    trait_h!(img, px0, px1, py0, axes_u8, ep_axe)
    trait_h!(img, px0, px1, py1, axes_u8, ep_axe)
    trait_v!(img, px0, py0, py1, axes_u8, ep_axe)

    if ligne_zero
        rz = rangee(0.0)
        for i in px0:(4 * T):px1                      # pointilles
            trait_h!(img, i, min(i + 2 * T, px1), rz, axes_u8, ep_axe)
        end
    end

    base = rangee(0.0)
    d = ep ÷ 2
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

    # age d'une colonne, en nombre de colonnes derriere la tete de lecture
    col(k) = mod(c_rev - 1 - k, n_tot) + 1
    alpha(k) = clamp(1.0 - k / n_tot, 0.0, 1.0)

    # deux passes du plus vieux au plus recent : les remplissages d'abord,
    # sinon le pinceau large d'une colonne pleine efface le front qui precede
    if n_dessin > 0 && c_rev > 0
        if remplir
            for k in (n_dessin - 1):-1:0
                c = col(k)
                x = px0 + c - 1
                rect!(img, x, x, min(base, rlo[c]), max(base, rlo[c]),
                      teinte(0.30 * alpha(k)))
            end
        end
        # pinceau carre : meme epaisseur sur les paliers et sur les fronts
        for k in (n_dessin - 1):-1:0
            c = col(k)
            x = px0 + c - 1
            rect!(img, x - d, x - d + ep - 1,
                  clamp(rlo[c] - d, py0, py1),
                  clamp(rhi[c] - d + ep - 1, py0, py1), teinte(alpha(k)))
        end
    end

    # tete de lecture
    if 0 < c_rev < length(donnees)
        trait_v!(img, px0 + c_rev - 1, py0, py1, axes_u8, max(2, ep ÷ 2))
    end

    ecris_texte!(img, titre, px0, haut + 2 * T, T, c_u8)
    ecris_texte!(img, etiq_haut,
                 px0 - largeur_texte(etiq_haut, T) - 2 * T, py0 - 3 * T, T, axes_u8)
    ecris_texte!(img, etiq_bas,
                 px0 - largeur_texte(etiq_bas, T) - 2 * T, py1 - 3 * T, T, axes_u8)
    return
end

"Peint le panneau de droite : x, y et puissance laser empiles verticalement."
function dessine_plots!(img, W_img::Int, H_img::Int, zones, traces,
                        c_rev::Int, n_dessin::Int, ep::Int,
                        t_a::Float64, t_b::Float64)
    fond_u8 = octets(FOND)
    axes_u8 = octets(COULEUR_AXES)
    ep_axe = max(2, ep ÷ 4)
    T = TAILLE_TEXTE
    rect!(img, W_img + 1, 2 * W_img, 1, H_img, fond_u8)
    trait_v!(img, W_img + 1, 1, H_img, axes_u8, ep_axe)

    gx, gy, g1, g2 = traces
    dessine_plot!(img, zones[1], gx, XLIM[1], XLIM[2], COULEUR_PLOT_X,
                  "x (mV)", @sprintf("%.0f", XLIM[2]), @sprintf("%.0f", XLIM[1]),
                  c_rev, n_dessin, ep, ep_axe, false, axes_u8, true)
    dessine_plot!(img, zones[2], gy, YLIM[1], YLIM[2], COULEUR_PLOT_Y,
                  "y (mV)", @sprintf("%.0f", YLIM[2]), @sprintf("%.0f", YLIM[1]),
                  c_rev, n_dessin, ep, ep_axe, false, axes_u8, true)
    dessine_plot!(img, zones[3], g1, -0.08, 1.08, COULEUR_PLOT_P1,
                  "puissance laser 1", "1", "0",
                  c_rev, n_dessin, ep, ep_axe, true, axes_u8, false)
    dessine_plot!(img, zones[4], g2, -0.08, 1.08, COULEUR_PLOT_P2,
                  "puissance laser 2", "1", "0",
                  c_rev, n_dessin, ep, ep_axe, true, axes_u8, false)

    # axe du temps sous la bande du bas
    px0, px1, _, py1, _ = zones[4]
    g = @sprintf("%.0f ms", t_a)
    d = @sprintf("%.0f ms", t_b)
    ecris_texte!(img, g, px0, py1 + 3 * T, T, axes_u8)
    ecris_texte!(img, d, px1 - largeur_texte(d, T), py1 + 3 * T, T, axes_u8)
    return
end

# -----------------------------------------------------------------------------
#  Simulation
# -----------------------------------------------------------------------------

function simule()
    isempty(POSITIONS) && error("POSITIONS est vide : aucune tornade a simuler.")
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

    sigma_px = max(0.6, SIGMA_FAISCEAU * ECHELLE_PX)
    rayon_px = RAYON_SPOT_SIGMA * sigma_px
    tau = TRAINEE_MS
    decroissance = tau > 0 ? exp(-dt_reel / tau) : 0.0

    n_sub = clamp(ceil(Int, SUR_ECHANTILLONNAGE * dt_reel / dwell), 1, 256)
    dt_sub = dt_reel / n_sub

    mkpath(DOSSIER_SORTIE)
    nom = isempty(NOM_FICHIER) ?
          "tornades_" * Dates.format(now(), "yyyy-mm-dd_HHhMMmSS") * ".mp4" :
          NOM_FICHIER
    chemin = joinpath(DOSSIER_SORTIE, nom)

    println("Tornades          : $n_tor position(s) x $N_REPETITIONS repetition(s)")
    println("Sequence          : laser 1, laser 2, delai, a chaque position")
    @printf("Temps reel        : %.0f ms (%.2f ms par point)\n", duree_reelle_ms, dwell)
    @printf("Video             : %d x %d px, %.0f fps, %.1f s (ralenti %.0fx)\n",
            W_tot, Hpx, FPS, n_frames / FPS, RALENTI)
    println("Images            : $n_frames ($n_sub sous-pas par image)")
    println("Sortie            : $chemin")
    if isodd(W_tot) || isodd(Hpx)
        println("Note              : dimension impaire, ffmpeg ajoute 1 px de bord noir.")
    end
    println()

    H1 = zeros(Float32, W, Hpx)                    # remanence du laser 1
    H2 = zeros(Float32, W, Hpx)                    # remanence du laser 2
    img = Array{UInt8}(undef, 3, W_tot, Hpx)
    txt = (round(UInt8, 255 * red(COULEUR_TEXTE)),
           round(UInt8, 255 * green(COULEUR_TEXTE)),
           round(UInt8, 255 * blue(COULEUR_TEXTE)))
    en_ms = duree_reelle_ms < 1000.0

    x_px(x) = (x - XLIM[1]) * ECHELLE_PX + 1.0
    y_px(y) = (YLIM[2] - y) * ECHELLE_PX + 1.0

    # graphiques : geometrie, echantillonnage et trace fixe si pas de fenetre
    zones = zones_plots(W, Hpx)
    n_col = zones[1][2] - zones[1][1] + 1
    ep_trace = EPAISSEUR_TRACE > 0 ? EPAISSEUR_TRACE : max(4, round(Int, W / 250))
    defilant = FENETRE_PLOTS_MS > 0
    duree_boucle_ms = n_tor * t_unite        # les plots montrent une boucle
    n_ech = clamp(ceil(Int, 4 * ((defilant ? FENETRE_PLOTS_MS : duree_boucle_ms) /
                                 n_col) / dwell), 2, 64)
    traces_fixes = (AFFICHER_PLOTS && !defilant) ?
        colonnes_traces(0.0, duree_boucle_ms, n_col, traj1, traj2, t_unite,
                        dwell, n_ech) : nothing

    taille = "$(W_tot)x$(Hpx)"
    vf = "pad=ceil(iw/2)*2:ceil(ih/2)*2:0:0:black"

    FFMPEG_jll.ffmpeg() do exe
        cmd = `$exe -hide_banner -loglevel error -y
               -f rawvideo -pix_fmt rgb24 -s $taille -r $FPS -i pipe:0
               -an -vf $vf -c:v libx264 -preset medium -crf $QUALITE_H264
               -pix_fmt yuv420p -movflags +faststart $chemin`

        open(cmd, "w") do flux
            for k in ProgressBar(1:n_frames)
                t1 = k * dt_reel                      # fin de l'image, temps reel

                # la remanence deja presente vieillit d'une image
                if decroissance == 0.0
                    fill!(H1, 0.0f0); fill!(H2, 0.0f0)
                else
                    H1 .*= Float32(decroissance); H2 .*= Float32(decroissance)
                end

                # exposition : sous-pas integres sur la duree de l'image
                for s in 1:n_sub
                    t = t1 - dt_reel + (s - 0.5) * dt_sub
                    x, y, laser = etat_faisceau(t, traj1, traj2, t_unite, dwell)
                    laser == 0 && continue
                    age = t1 - t
                    poids = (dt_sub / dwell) * (tau > 0 ? exp(-age / tau) : 1.0)
                    depose_spot!(laser == 1 ? H1 : H2,
                                 x_px(x), y_px(y), sigma_px, poids, rayon_px)
                end

                colorise!(img, H1, H2, COULEUR_LASER_1, COULEUR_LASER_2, FOND)

                if AFFICHER_TEMPS
                    etiq = en_ms ? @sprintf("t = %6.1f ms", t1) :
                                   @sprintf("t = %6.3f s", t1 / 1000)
                    marge = max(12, round(Int, 0.025 * W))
                    ecris_texte!(img, etiq, marge, marge, TAILLE_TEXTE, txt)
                end
                AFFICHER_ECHELLE && barre_echelle!(img, W, ECHELLE_PX, txt)

                if AFFICHER_PLOTS
                    if defilant
                        ta, tb = t1 - FENETRE_PLOTS_MS, t1
                        tr = colonnes_traces(ta, tb, n_col, traj1, traj2,
                                             t_unite, dwell, n_ech)
                        dessine_plots!(img, W, Hpx, zones, tr, n_col, n_col,
                                       ep_trace, ta, tb)
                    else
                        # la tete de lecture repart de zero a chaque repetition
                        r = t1 - floor(t1 / duree_boucle_ms) * duree_boucle_ms
                        premier = t1 < duree_boucle_ms
                        c_rev = clamp(round(Int, r / duree_boucle_ms * n_col),
                                      0, n_col)
                        # pile sur la couture le reste tombe a ~0 : la tete est
                        # en fin de parcours, pas avant le debut
                        !premier && c_rev == 0 && (c_rev = n_col)
                        # apres le premier tour, le trace de la boucle
                        # precedente reste visible et finit de s'effacer
                        n_dessin = premier ? c_rev : n_col
                        dessine_plots!(img, W, Hpx, zones, traces_fixes, c_rev,
                                       n_dessin, ep_trace, 0.0, duree_boucle_ms)
                    end
                end

                write(flux, img)
            end
        end
    end

    println("\nVideo ecrite : $chemin")
    return chemin
end

simule()
