using Plots
using Random

# ============================== RÉGLAGES ==============================

baseline   = 0.0      # valeur au repos (avant le step)
steady     = 1.0      # valeur asymptotique pendant le step (à t = t_on)
tau_on     = 0.3      # constante de temps de la montée (réponse au step)
tau_off    = 1.0      # constante de temps de la descente (relâche)

t_start = 0.0          # début de la simulation
t_on    = 1.0          # instant où le step est appliqué
t_off   = 4.0           # instant où le step est relâché
t_end   = 8.0           # fin de la simulation

dt = 0.001

noise_std = 0.25        # écart-type du bruit gaussien ajouté au signal

# Couleurs du fond et des axes
background_color = :black
axis_color = :white

signal_color = :dodgerblue
clean_color = :white

step_color = :gold
step_alpha = 0.2

# ========================================================================

t = t_start:dt:t_end

y = similar(collect(t))
for (i, ti) in enumerate(t)
    if ti < t_on
        y[i] = baseline
    elseif ti < t_off
        y[i] = baseline + (steady - baseline) * (1 - exp(-(ti - t_on) / tau_on))
    else
        y_at_off = baseline + (steady - baseline) * (1 - exp(-(t_off - t_on) / tau_on))
        y[i] = baseline + (y_at_off - baseline) * exp(-(ti - t_off) / tau_off)
    end
end

y_noisy = y .+ noise_std .* randn(length(y))

plt = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

y_min, y_max = extrema(y_noisy)
step_band = Shape([t_on, t_off, t_off, t_on], [y_min, y_min, y_max, y_max])
plot!(plt, step_band, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

plot!(plt, t, y_noisy, color = signal_color, linewidth = 1, alpha = 0.6)
plot!(plt, t, y, color = clean_color, linewidth = 2)

display(plt)

# ============================== SÉRIE DE STEPS ==============================

step_levels    = [0, 1, 0, 2, 0, 3, 0, 4, 0]  # niveaux successifs de l'entrée
step_duration  = 2.0                          # durée (temps) de chaque niveau
tau_on         = 0.2      # constante de temps de la montée (réponse au step)
tau_off        = 0.3      # constante de temps de la descente (relâche)

t2 = 0:dt:(step_duration * length(step_levels))

# ---- calcul de la réponse du premier ordre à la série de steps ----
y2 = similar(collect(t2))
y_current = Float64(step_levels[1])
for (i, ti) in enumerate(t2)
    seg = min(floor(Int, ti / step_duration) + 1, length(step_levels))
    level = Float64(step_levels[seg])
    t_seg_start = (seg - 1) * step_duration
    if seg == 1
        y2[i] = level
    else
        y_seg_start = Float64(step_levels[seg - 1])
        tau = level >= y_seg_start ? tau_on : tau_off
        y2[i] = level + (y_seg_start - level) * exp(-(ti - t_seg_start) / tau)
    end
end

y2_noisy = y2 .+ noise_std .* randn(length(y2))

plt2 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

y2_min, y2_max = extrema(y2_noisy)
for (seg, level) in enumerate(step_levels)
    level == 0 && continue
    t_seg_start = (seg - 1) * step_duration
    t_seg_end = seg * step_duration
    band = Shape([t_seg_start, t_seg_end, t_seg_end, t_seg_start], [y2_min, y2_min, y2_max, y2_max])
    plot!(plt2, band, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)
end

plot!(plt2, t2, y2_noisy, color = signal_color, linewidth = 1, alpha = 0.6)
plot!(plt2, t2, y2, color = clean_color, linewidth = 2)

display(plt2)

# ============================== SETPOINTS BASAUX DIFFÉRENTS ==============================

basal_levels = [0.0, -0.25, -0.5, -0.75]   # setpoint basal de chaque courbe
curve_colors = [:dodgerblue, :seagreen, :orchid, :orangered]

tau_on     = 0.3      # constante de temps de la montée (réponse au step)
tau_off    = 1.0      # constante de temps de la descente (relâche)
noise_std = 0.1        # écart-type du bruit gaussien ajouté au signal

# ========================================================================

y3s = [similar(collect(t)) for _ in basal_levels]
for (k, basal) in enumerate(basal_levels)
    y3 = y3s[k]
    for (i, ti) in enumerate(t)
        if ti < t_on
            y3[i] = basal
        elseif ti < t_off
            y3[i] = basal + (steady - basal) * (1 - exp(-(ti - t_on) / tau_on))
        else
            y_at_off = basal + (steady - basal) * (1 - exp(-(t_off - t_on) / tau_on))
            y3[i] = basal + (y_at_off - basal) * exp(-(ti - t_off) / tau_off)
        end
    end
end

y3s_noisy = [y3 .+ noise_std .* randn(length(y3)) for y3 in y3s]

plt3 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

y3_min = minimum(minimum.(y3s_noisy))
y3_max = maximum(maximum.(y3s_noisy))
step_band3 = Shape([t_on, t_off, t_off, t_on], [y3_min, y3_min, y3_max, y3_max])
plot!(plt3, step_band3, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

for (k, color) in enumerate(curve_colors)
    plot!(plt3, t, y3s_noisy[k], color = color, linewidth = 1, alpha = 0.6)
end
for (k, color) in enumerate(curve_colors)
    plot!(plt3, t, y3s[k], color = color, linewidth = 2)
end

display(plt3)

# ============================== VALEURS DE STEP DIFFÉRENTES ==============================

step_values  = [0.1, 1.0]   # amplitude du step pour chaque courbe (même basal)
curve_colors = [:dodgerblue, :orangered]

# ========================================================================

y4s = [similar(collect(t)) for _ in step_values]
for (k, step_value) in enumerate(step_values)
    y4 = y4s[k]
    for (i, ti) in enumerate(t)
        if ti < t_on
            y4[i] = baseline
        elseif ti < t_off
            y4[i] = baseline + (step_value - baseline) * (1 - exp(-(ti - t_on) / tau_on))
        else
            y_at_off = baseline + (step_value - baseline) * (1 - exp(-(t_off - t_on) / tau_on))
            y4[i] = baseline + (y_at_off - baseline) * exp(-(ti - t_off) / tau_off)
        end
    end
end

y4s_noisy = [y4 .+ noise_std .* randn(length(y4)) for y4 in y4s]

plt4 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

y4_min = minimum(minimum.(y4s_noisy))
y4_max = maximum(maximum.(y4s_noisy))
step_band4 = Shape([t_on, t_off, t_off, t_on], [y4_min, y4_min, y4_max, y4_max])
plot!(plt4, step_band4, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

for (k, color) in enumerate(curve_colors)
    plot!(plt4, t, y4s_noisy[k], color = color, linewidth = 1, alpha = 0.6)
end
for (k, color) in enumerate(curve_colors)
    plot!(plt4, t, y4s[k], color = color, linewidth = 2)
end

display(plt4)

# ============================== SETPOINTS ET STEPS DIFFÉRENTS ==============================

basal_levels5 = [0.0, 0.1, 0.2, 0.3]   # setpoint basal de chaque courbe
step_values5  = [1.0, 1.5, 2.0, 2.5]      # amplitude du step pour chaque courbe
curve_colors = [:dodgerblue, :seagreen, :orchid, :orangered]

# ========================================================================

y5s = [similar(collect(t)) for _ in basal_levels5]
for (k, (basal, step_value)) in enumerate(zip(basal_levels5, step_values5))
    y5 = y5s[k]
    for (i, ti) in enumerate(t)
        if ti < t_on
            y5[i] = basal
        elseif ti < t_off
            y5[i] = basal + (step_value - basal) * (1 - exp(-(ti - t_on) / tau_on))
        else
            y_at_off = basal + (step_value - basal) * (1 - exp(-(t_off - t_on) / tau_on))
            y5[i] = basal + (y_at_off - basal) * exp(-(ti - t_off) / tau_off)
        end
    end
end

y5s_noisy = [y5 .+ noise_std .* randn(length(y5)) for y5 in y5s]

plt5 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

y5_min = minimum(minimum.(y5s_noisy))
y5_max = maximum(maximum.(y5s_noisy))
step_band5 = Shape([t_on, t_off, t_off, t_on], [y5_min, y5_min, y5_max, y5_max])
plot!(plt5, step_band5, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

for (k, color) in enumerate(curve_colors)
    plot!(plt5, t, y5s_noisy[k], color = color, linewidth = 1, alpha = 0.6)
end
for (k, color) in enumerate(curve_colors)
    plot!(plt5, t, y5s[k], color = color, linewidth = 2)
end

display(plt5)

# ============================== CHLORURE (6) ET CALCIUM (7) ==============================
#
# Plot 6 : [Cl⁻]i intracellulaire — reprend exactement la courbe du 1er graphique.
# Plot 7 : [Ca²⁺]i en réponse à cette variation de [Cl⁻]i, en régime
#          « spike-triggered » : chaque potentiel d'action déclenche un
#          transitoire calcique local rapide (montée < 1 ms). Le chlorure ne
#          dessine plus directement la courbe : il module la FRÉQUENCE des PA
#          (voir plus bas).
#
# ---- Justification biologique (mécanisme direct du Cl⁻, sans GABA) ----
#
# Le Cl⁻ est un anion : plus sa concentration intracellulaire est élevée, plus
# la cellule porte de charge négative mobilisable pour hyperpolariser/shunter
# la membrane dès qu'une conductance Cl⁻ est active. Cet effet ÉLOIGNE le
# potentiel de membrane du seuil de décharge :
#
#   [Cl⁻]i ↑  →  hyperpolarisation / shunt de la membrane
#            →  seuil de PA plus difficile à atteindre
#            →  MOINS de potentiels d'action.
#
# On n'invoque plus ici le switch dépolarisant du GABA_A (E_GABA, efflux,
# etc.) : le Cl⁻ agit dans son rôle inhibiteur direct et classique. Chaque PA
# qui parvient tout de même à se déclencher reste un évènement tout-ou-rien :
# il produit un transitoire calcique local rapide de taille FIXE (montée
# < 1 ms) — seule la FRÉQUENCE des PA dépend du chlorure, pas l'amplitude.

# ---- ÉCHELLE DE TEMPS ----
# Les constantes du chlorure sont fixées (plot 6 = plot 1) : tau_on = 0.3,
# tau_off = 1.0. On les calibre sur la récupération du Cl⁻ par KCC2, mesurée à
# τ ≈ 6.5-7.2 s dans le néocortex de rat (≈ 18-19 s chez l'humain).
#   →  1 unité de temps ≈ 7 s
# TOUS les paramètres ci-dessous sont vérifiés contre cette échelle.

SEC_PER_UNIT = 7.0    # 1 unité de temps ≈ 7 s

cl_tau_on  = 0.3      # identiques au plot 1 (définis ici en dur, indépendants
cl_tau_off = 1.0      # de tau_on/tau_off qui sont réutilisés par les plots
                      # précédents et ne pointent plus vers 0.3/1.0 à ce stade)

Ca_base7 = 60.0       # nM, [Ca²⁺]i de repos
tau_ca7  = 0.15       # s, tamponnage + extrusion entre les PA. Le transitoire
                      # reste LOCAL (bouton/dendrite), pas somatique : la plage
                      # mesurée va de 20-40 ms (épines) à ~80 ms (dendrites
                      # proximales, 35-37 °C) jusqu'à ~260-290 ms (composante
                      # lente des boutons synaptiques) — 150 ms se situe dans
                      # cette dernière plage, plus lent qu'avant tout en
                      # restant bien plus rapide que la constante somatique
                      # globale (0.37-1.7 s) utilisée ailleurs pour un signal
                      # en volume.

# ========================================================================

# ---- Plot 6 : chlorure (identique au 1er graphique) ----
y_cl = y
y_cl_noisy = y_noisy

plt6 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

cl_min, cl_max = extrema(y_cl_noisy)
band6 = Shape([t_on, t_off, t_off, t_on], [cl_min, cl_min, cl_max, cl_max])
plot!(plt6, band6, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

plot!(plt6, t, y_cl_noisy, color = signal_color, linewidth = 1, alpha = 0.9)
plot!(plt6, t, y_cl, color = clean_color, linewidth = 2)

display(plt6)

# ---- Plot 7 : calcium « spike-triggered » en réponse au chlorure ----
#
# Les PA sont tirés AU HASARD dans le temps (processus de Poisson inhomogène
# avec période réfractaire), donc répartis sur toute la fenêtre mais jamais à
# intervalles réguliers.

t7 = collect(t_start:dt:t_end)   # même fenêtre que le chlorure des plots 1/6
dts7 = dt * SEC_PER_UNIT

# --- chlorure : même génération qu'au plot 1, normalisée 0 (repos) -> 1 (charge) ---
cl_drive7 = similar(t7)
for (i, ti) in enumerate(t7)
    if ti < t_on
        cl_drive7[i] = 0.0
    elseif ti < t_off
        cl_drive7[i] = 1 - exp(-(ti - t_on) / cl_tau_on)
    else
        y_at = 1 - exp(-(t_off - t_on) / cl_tau_on)
        cl_drive7[i] = y_at * exp(-(ti - t_off) / cl_tau_off)
    end
end

Random.seed!(12)

rate_rest7 = 2.2     # Hz, fréquence de décharge spontanée au repos ([Cl⁻]i bas)
rate_min7  = 0.2     # Hz, fréquence résiduelle à la charge maximale de [Cl⁻]i
                     # (fortement réduite : effet inhibiteur direct du Cl⁻ —
                     # même contraste relatif ~11x qu'avant, juste plus de PA
                     # en valeur absolue des deux côtés)
amp7       = 55.0    # nM, saut de Ca²⁺ par PA — FIXE : un PA est tout-ou-rien,
                     # sa taille ne dépend pas du chlorure ici
refrac7    = 0.01    # s, période réfractaire
ca_noise7  = 3.5     # nM, bruit de mesure

tau_ca_u7 = tau_ca7 / SEC_PER_UNIT   # extrusion, en unités de temps du tracé

ca7 = similar(t7)
spikes7 = Float64[]
let c = Ca_base7, last_spk = -Inf
    for i in eachindex(t7)
        # relation INVERSE : plus de Cl⁻ -> moins de PA
        rate = rate_rest7 - (rate_rest7 - rate_min7) * cl_drive7[i]

        # extrusion : relaxation vers la ligne de base
        c += dt * (Ca_base7 - c) / tau_ca_u7

        # tirage de Poisson : un PA -> saut instantané de Ca²⁺
        if (t7[i] - last_spk) * SEC_PER_UNIT > refrac7 && rand() < rate * dts7
            c += amp7
            push!(spikes7, t7[i])
            last_spk = t7[i]
        end

        ca7[i] = c
    end
end

ca7_noisy = ca7 .+ ca_noise7 .* randn(length(ca7))

plt7 = plot(
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

ca7_min, ca7_max = extrema(ca7_noisy)
band7 = Shape([t_on, t_off, t_off, t_on], [ca7_min, ca7_min, ca7_max, ca7_max])
plot!(plt7, band7, fillcolor = step_color, fillalpha = step_alpha, linealpha = 0)

hline!(plt7, [Ca_base7], linestyle = :dot, color = :gray, linealpha = 0.7)

plot!(plt7, t7, ca7_noisy, color = :orangered, linewidth = 1, alpha = 0.5)
# plot!(plt7, t7, ca7, color = clean_color, linewidth = 1.2)

display(plt7)
