using LibSerialPort
using Plots
using Random
using LinearAlgebra
using Statistics
using QuadGK
using SpecialFunctions
using ProgressBars

try
    global serial = LibSerialPort.open("/dev/tty.usbmodem11102", 115200)
catch e
    println("Error opening serial port: $e")
end

function send_command(serial, command)
    write(serial, command)
    sleep(0.000001)
    response = strip(read(serial, String))
    return response
end

function get_response(serial)
    response = strip(read(serial, String))
    return response
end

function send_points_as_ab(serial, pts::Vector{Tuple{Int, Int}}, f::Real)
    response = send_command(serial, "C Z\n")

    if occursin("ERR", response)
        println("Error ($response) in command C Z")
    end

    for (idx, (x, y)) in ProgressBar(enumerate(pts))
        response = send_command(serial, "A 0 AB 1 $idx $x\n")
        if occursin("ERR", response)
            println("Error ($response) in command A 0 AB 1 $idx $x")
        end
        response = send_command(serial, "A 0 AB 2 $idx $y\n")
        if occursin("ERR", response)
            println("Error ($response) in command A 0 AB 2 $idx $y")
        end
    end
end

function gaussian_power_inside_disk(ρ::Real, R::Real, σ::Real)
    integrand(r) = (r / σ^2) * exp(-(r^2 + ρ^2) / (2σ^2)) * besseli(0, (r * ρ) / σ^2)
    val, _ = quadgk(integrand, 0.0, R)
    return val
end

function compensated_radial_cdf(R::Real, σ::Real; M::Int = 500)
    ρgrid = range(0.0, R, length = M)
    w = zeros(Float64, M)

    for i in 1:M
        f = gaussian_power_inside_disk(ρgrid[i], R, σ)
        w[i] = 1.0 / max(f, 1e-12)
    end

    cdf = zeros(Float64, M)

    for i in 2:M
        dρ = ρgrid[i] - ρgrid[i - 1]
        cdf[i] = cdf[i - 1] + 0.5 * (
            ρgrid[i] * w[i] + ρgrid[i - 1] * w[i - 1]
        ) * dρ
    end

    cdf ./= cdf[end]

    return collect(ρgrid), cdf
end

function invert_cdf(ρgrid::Vector{Float64}, cdf::Vector{Float64}, u::Real)
    i = searchsortedfirst(cdf, u)

    if i <= 1
        return ρgrid[1]
    elseif i > length(cdf)
        return ρgrid[end]
    end

    u1, u2 = cdf[i - 1], cdf[i]
    ρ1, ρ2 = ρgrid[i - 1], ρgrid[i]

    t = (u - u1) / max(u2 - u1, 1e-12)

    return ρ1 + t * (ρ2 - ρ1)
end

function spiral_points(N::Int, R::Real, σ::Real;
                                   center::Tuple{Real,Real} = (0.0, 0.0),
                                   turns::Real = 8.0,
                                   M::Int = 500)

    x0, y0 = float(center[1]), float(center[2])
    ρgrid, cdf = compensated_radial_cdf(R, σ; M = M)

    points = Vector{Tuple{Int64,Int64}}(undef, N)

    for k in 1:N
        u = (k - 0.5) / N
        ρ = invert_cdf(ρgrid, cdf, u)
        θ = π * turns * u

        x = x0 + ρ * cos(θ)
        y = y0 + ρ * sin(θ)

        points[k] = (round(Int64, x), round(Int64, y))
    end

    return points
end



dist2(p, q) = hypot(float(q[1]) - float(p[1]), float(q[2]) - float(p[2]))

function path_length_cycle(points::AbstractVector{<:Tuple{<:Real,<:Real}})
    n = length(points)
    n <= 1 && return 0.0

    s = 0.0
    for i in 1:n-1
        s += dist2(points[i], points[i+1])
    end
    s += dist2(points[end], points[1])
    return s
end

function nearest_neighbor_cycle(points::AbstractVector{<:Tuple{<:Real,<:Real}})
    n = length(points)
    n <= 1 && return collect(points)

    used = falses(n)
    order = Vector{Int}(undef, n)

    current = 1
    order[1] = current
    used[current] = true

    for k in 2:n
        best_j = 0
        best_d = Inf
        p = points[current]

        for j in 1:n
            if !used[j]
                d = dist2(p, points[j])
                if d < best_d
                    best_d = d
                    best_j = j
                end
            end
        end

        order[k] = best_j
        used[best_j] = true
        current = best_j
    end

    return [points[i] for i in order]
end

function two_opt_cycle!(tour::Vector{Tuple{Float64,Float64}})
    n = length(tour)
    n <= 3 && return tour

    improved = true
    while improved
        improved = false

        for i in 2:n-2
            for k in i+1:n-1
                A = tour[i-1]
                B = tour[i]
                C = tour[k]
                D = tour[k+1]

                old = dist2(A, B) + dist2(C, D)
                new = dist2(A, C) + dist2(B, D)

                if new + 1e-12 < old
                    reverse!(tour, i, k)
                    improved = true
                end
            end
        end
    end

    return tour
end

function heuristic_tour(points::AbstractVector{<:Tuple{<:Real,<:Real}})
    tour = [ (float(p[1]), float(p[2])) for p in nearest_neighbor_cycle(points) ]
    two_opt_cycle!(tour)
    return tour
end

function show_progress(current::Int, total::Int; width::Int = 40, prefix::String = "Progress")
    ratio = total == 0 ? 1.0 : current / total
    filled = clamp(round(Int, ratio * width), 0, width)
    bar = repeat("█", filled) * repeat(" ", width - filled)
    pct = round(100 * ratio; digits = 1)

    print("\r$prefix [$bar] $(lpad(string(pct), 5))% ($current/$total)")
    flush(stdout)

    if current >= total
        println()
    end
end

function optimize_centers(points::AbstractVector{<:Tuple{<:Real,<:Real}})
    n = length(points)
    n <= 1 && return collect(points)

    pts = [(float(p[1]), float(p[2])) for p in points]

    # Borne supérieure initiale donnée par une heuristique
    best_guess = heuristic_tour(pts)
    upper_bound = path_length_cycle(best_guess)

    # Matrice des distances
    d = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        d[i, j] = dist2(pts[i], pts[j])
    end

    # On fixe le point 1 pour casser les symétries
    N = n - 1
    total_masks = 1 << N

    # dp[mask+1, j] = coût minimal pour partir de 1,
    # visiter exactement mask (sur les points 2..n),
    # et finir en j
    dp = fill(Inf, total_masks, n)
    parent = fill(UInt16(0), total_masks, n)

    # Cas de base : un seul point visité
    for j in 2:n
        mask = 1 << (j - 2)
        dp[mask + 1, j] = d[1, j]
        parent[mask + 1, j] = UInt16(1)
    end

    # Petite optimisation : on suit la progression sur les sous-ensembles
    step = max(1, total_masks ÷ 200)

    # Remplissage du DP
    for mask in ProgressBar(0:total_masks-1)

        # Pour chaque point j présent dans mask
        for j in 2:n
            bitj = 1 << (j - 2)
            if (mask & bitj) == 0
                continue
            end

            cur = dp[mask + 1, j]
            if !isfinite(cur) || cur >= upper_bound
                continue
            end

            # Extension vers un point k non encore visité
            remaining = (~mask) & (total_masks - 1)
            while remaining != 0
                lb = remaining & -remaining
                kbit = trailing_zeros(lb)
                k = kbit + 2
                newmask = mask | lb
                newcost = cur + d[j, k]

                if newcost < dp[newmask + 1, k] && newcost < upper_bound
                    dp[newmask + 1, k] = newcost
                    parent[newmask + 1, k] = UInt16(j)
                end

                remaining -= lb
            end
        end
    end

    # Fermeture du cycle
    fullmask = total_masks - 1
    best_cost = Inf
    best_last = 0

    for j in 2:n
        c = dp[fullmask + 1, j] + d[j, 1]
        if c < best_cost
            best_cost = c
            best_last = j
        end
    end

    # Reconstruction du chemin optimal
    order = Vector{Int}(undef, n)
    order[1] = 1

    mask = fullmask
    last = best_last

    for pos in n:-1:2
        order[pos] = last
        prev = Int(parent[mask + 1, last])
        mask = mask & ~(1 << (last - 2))
        last = prev
    end

    return [points[i] for i in order]
end




i = 100
t = 90
dt = 10
j = round(Int64, dt/t * i)

n = 10
r = 100
σ = 10

x0 = 0
y0 = 0

m = 10 # maximum 20 pour le moment
f = 1000/((t + dt)*m)

centers = [(rand(-1000:1000), rand(-1000:1000)) for _ in 1:m]

println("Détermination du path optimal...")
centers = optimize_centers(centers)
println("Path optimal déterminé.")

# pts = spiral_points_uniform(100, 500, center=(0,0), turns=4)

pts = Vector{Tuple{Int64, Int64}}([])

# append!(pts, repeat([(0, 0)], j))
# append!(pts, spiral_points(i, r, σ, center=(x0, y0), turns=n))

first = true

for center in centers
    if first 
        global first = false
    else
        append!(pts, repeat([center], j))
    end

    append!(pts, spiral_points(i, r, σ, center=center, turns=n))
end

append!(pts, repeat([centers[1]], j))


# append!(pts, spiral_points(i, r, σ, center=(0, 1000), turns=n))
# append!(pts, spiral_points(i, r, σ, center=(-1000, 1000), turns=n))
# append!(pts, spiral_points(i, r, σ, center=(-1000, 0), turns=n))
# pts = spiral_points_for_gaussian(500, 0, 0; npts=100, turns=4, spot_sigma_mv=100.0, step_factor=0.nnonon5, dense=10000)

println("Envoi des commandes au buffer...")
send_points_as_ab(serial, pts, 1)
println("Commandes envoyées au buffer.")

response = send_command(serial, "S 1 D 1 R 1\n")
response *= send_command(serial, "A 1 AA $f $(length(pts))\n")
response *= send_command(serial, "A 0 DD 2 $(f*m) $(100-round(Int64, 100*dt/(dt + t))) 1 $m 10\n")

if occursin("ERR", response)
    println("Error ($response) in command A 0 AA $f $(length(pts))")
end

LibSerialPort.close(serial)

x_pts = [p[1] for p in pts]
y_pts = [p[2] for p in pts]
scatter(x_pts, y_pts, title="Tornado Points", xlabel="X (mV)", ylabel="Y (mV)", legend=false, aspect_ratio=:equal)
