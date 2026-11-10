# regulateur.jl — PI avec anti-emballement, sur l'estimation du filtre.
#
#   u_k = sat_[0,1]( Kp·e_k + I_k ),  e_k = C_consigne − Ĉ_k,  I_{k+1} = I_k + Ki·e_k·Δt
#
# Anti-emballement par intégration conditionnelle : l'intégrateur est gelé quand la commande est
# saturée et que l'écart pousse dans le sens de la saturation (u = 0 avec e < 0 : la consigne
# est sous le chlore, seul KCC2 peut le faire descendre). Pas de terme dérivé.
#
# Trois variantes, selon ce que le régulateur fait de P₁₁ (corrigé par c = P₁₁/err² mesuré) :
#   :A aveugle ; :B zone morte, ne corrige que si |e| > √(P₁₁/c) ;
#   :C gain pondéré Kp·σ₀²/(P₁₁/c + σ₀²), plus mou quand le filtre doute.

mutable struct PI
    Kp::Float64
    Ki::Float64
    I::Float64
    variante::Symbol
    c::Float64          # calibration du filtre : P₁₁ annoncé ÷ erreur quadratique réelle
    σ0²::Float64        # variante C
end

PI(Kp, Ki; variante = :A, c = 1.0, σ0² = 1.0) = PI(Kp, Ki, 0.0, variante, c, σ0²)

"""La commande du pas, à partir de l'écart e et de P₁₁ ; met l'intégrateur à jour."""
function commande!(r::PI, e, P11, dt)
    Kp, Ki = r.Kp, r.Ki
    if r.variante === :B
        abs(e) > sqrt(P11 / r.c) || (e = 0.0)
    elseif r.variante === :C
        w = r.σ0² / (P11 / r.c + r.σ0²)
        Kp, Ki = Kp * w, Ki * w
    end
    v = Kp * e + r.I
    u = clamp(v, 0.0, 1.0)
    libre = (0 < v < 1) || (v >= 1 && e < 0) || (v <= 0 && e > 0)
    libre && (r.I += Ki * e * dt)
    return u
end
