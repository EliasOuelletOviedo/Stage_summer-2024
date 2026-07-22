using LibSerialPort
using Plots
using Random
using LinearAlgebra
using Statistics
using QuadGK
using SpecialFunctions
using ProgressBars
using Images
using Combinatorics

include("ImageJROI.jl")
using .ImageJROI

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

function send_points(serial, pts::Vector{Tuple{Int, Int}})
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

# Distance depuis (cx, cy) jusqu'à l'intersection entre le rayon d'angle θ
# et le contour polygonal (xs, ys). Retourne 0.0 si aucune intersection trouvée.
function polygon_ray_distance(θ::Real, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real},
                                          cx::Real, cy::Real)
    dx, dy = cos(θ), sin(θ)
    n = length(xs)
    best = Inf

    for i in 1:n
        j = i == n ? 1 : i + 1
        ax, ay = xs[i] - cx, ys[i] - cy
        ex, ey = xs[j] - xs[i], ys[j] - ys[i]

        det = ex * dy - ey * dx
        abs(det) < 1e-12 && continue

        t = (ex * ay - ey * ax) / det
        s = (dx * ay - dy * ax) / det

        if t > 1e-9 && 0.0 <= s <= 1.0 && t < best
            best = t
        end
    end

    return isfinite(best) ? best : 0.0
end

# Variante de spiral_points qui épouse la forme du ROI (xs, ys) plutôt qu'un
# disque : la densité radiale reste calquée sur compensated_radial_cdf (calculée
# pour un rayon moyen R), mais chaque point est ensuite mis à l'échelle vers la
# distance réelle du contour dans la direction θ.
function shape_spiral_points(N::Int, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, σ::Real;
                                         center::Tuple{Real,Real} = (0.0, 0.0),
                                         turns::Real = 8.0,
                                         M::Int = 500)

    x0, y0 = float(center[1]), float(center[2])

    R = mean(hypot(xs[k] - x0, ys[k] - y0) for k in eachindex(xs))
    ρgrid, cdf = compensated_radial_cdf(R, σ; M = M)

    points = Vector{Tuple{Int64,Int64}}(undef, N)

    for k in 1:N
        u = (k - 0.5) / N
        ρ = invert_cdf(ρgrid, cdf, u)
        θ = π * turns * u

        R_θ = polygon_ray_distance(θ, xs, ys, x0, y0)
        R_θ = R_θ > 0 ? R_θ : R
        ρ_scaled = ρ * (R_θ / R)

        x = x0 + ρ_scaled * cos(θ)
        y = y0 + ρ_scaled * sin(θ)

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

function two_opt_cycle!(tour::Vector{Tuple{Float32,Float32}})
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



function circle_from_2(a, b)
    center = ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)
    radius = hypot(a[1] - center[1], a[2] - center[2])
    return center, radius
end

function circle_from_3(a, b, c)
    A = [
        b[1]-a[1]  b[2]-a[2];
        c[1]-a[1]  c[2]-a[2]
    ]

    if abs(det(A)) < 1e-12
        return nothing
    end

    B = [
        ((b[1]^2 - a[1]^2) + (b[2]^2 - a[2]^2))/2;
        ((c[1]^2 - a[1]^2) + (c[2]^2 - a[2]^2))/2
    ]

    sol = A \ B
    center = (sol[1], sol[2])
    radius = hypot(a[1] - center[1], a[2] - center[2])

    return center, radius
end

function contains_all(center, radius, pts; eps=1e-9)
    for p in eachrow(pts)
        if hypot(p[1] - center[1], p[2] - center[2]) > radius + eps
            return false
        end
    end
    return true
end

function minimum_enclosing_circle(x_coords, y_coords)
    best_center = (0.0, 0.0)
    best_radius = Inf
    n = length(x_coords)

    for (i, j) in combinations(1:n, 2)
        c, r = circle_from_2((x_coords[i], y_coords[i]), (x_coords[j], y_coords[j]))
        if r < best_radius && contains_all(c, r, [x_coords y_coords])
            best_center = c
            best_radius = r
        end
    end

    for (i, j, k) in combinations(1:n, 3)
        result = circle_from_3((x_coords[i], y_coords[i]), (x_coords[j], y_coords[j]), (x_coords[k], y_coords[k]))
        result === nothing && continue

        c, r = result
        if r < best_radius && contains_all(c, r, [x_coords y_coords])
            best_center = c
            best_radius = r
        end
    end

    return best_center, best_radius
end

function centroid_circle(x_coords, y_coords)
    cx = mean(x_coords)
    cy = mean(y_coords)

    center = (cx, cy)

    radius = mean(
        hypot(p[1] - cx, p[2] - cy)
        for p in eachrow([x_coords y_coords])
    )

    return center, radius
end




filepath = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Test galvo/fov-001-C2-T005.tif"
img = load(filepath)
ny, nx = size(img)
plot(img, axis=true, ticks=true, legend=false, background_color = RGBA(18.0/255.0, 20.0/255.0, 20.0/255.0))

roi_folderpath = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Test galvo/RoiSet"
roi_paths = readdir(roi_folderpath; join=true)

zip = read_roi_zip("/Users/eliasouellet-oviedo/Documents/Maîtrise/Test galvo/RoiSet.zip")

i = 100
t = 950
dt = 50
j = round(Int64, dt/t * i)

n = 10
r = 100
σ = 10

x0 = 0
y0 = 0

m = length(roi_paths)
f = 1000/((t + dt)*m)

v_min_x, v_max_x = -1000, 1000
v_min_y, v_max_y = -500, 500

ROIs = []
pts = Vector{Tuple{Int64, Int64}}([])
buffer = Vector{Tuple{Int64, Int64}}([])
first = true

println("Extraction des ROI...")
for (name, roi) in ProgressBar(zip)
    x_coords = roi.x_coordinates
    y_coords = roi.y_coordinates
    
    plot!(x_coords, y_coords, color=:blue)

    # center, radius = minimum_enclosing_circle(x_coords, y_coords)
    center, radius = centroid_circle(x_coords, y_coords)

    push!(ROIs, [center, radius, x_coords, y_coords])
end
println("Extraction des ROI terminée.")

println("Détermination du path optimal...")
centers = [ROI[1] for ROI in ROIs]
centers = optimize_centers(centers)
ROIs = [ROIs[findfirst(x -> x[1] == center, ROIs)] for center in centers]
println("Path optimal déterminé.")

for ROI in ProgressBar(ROIs)
    center = round.(Int64, ROI[1])
    x_coords = ROI[3]
    y_coords = ROI[4]

    center_v = (-round(Int64, v_min_x + (center[1] - 1)*(v_max_x - v_min_x)/(nx - 1)),
                -round(Int64, v_min_y + (center[2] - 1)*(v_max_y - v_min_y)/(ny - 1)))

    x_coords_v = [-(v_min_x + (x - 1)*(v_max_x - v_min_x)/(nx - 1)) for x in x_coords]
    y_coords_v = [-(v_min_y + (y - 1)*(v_max_y - v_min_y)/(ny - 1)) for y in y_coords]

    if first
        global first_center = center_v
        global first = false
    else
        append!(buffer, repeat([center_v], j))
    end

    append!(pts, shape_spiral_points(i, x_coords, y_coords, σ, center=center, turns=n))
    append!(buffer, shape_spiral_points(i, x_coords_v, y_coords_v, σ, center=center_v, turns=n))
end

append!(buffer, repeat([first_center], j))

println(buffer)

println("Envoi des commandes au buffer...")
send_points(serial, buffer)
println("Commandes envoyées au buffer.")

response = send_command(serial, "S 1 D 1 R 1\n")
response *= send_command(serial, "A 1 AA $f $(length(buffer))\n")
response *= send_command(serial, "A 1 DD 2 $(f*m) $(100-round(Int64, 100*dt/(dt + t))) 1 $m 10\n")
response *= send_command(serial, "A 0 DO 3 1\n")

if occursin("ERR", response)
    println("Error ($response)")
end

LibSerialPort.close(serial)



x_pts = [p[1] for p in pts]
y_pts = [p[2] for p in pts]
scatter!(x_pts, y_pts, alpha = 0.1, color = :red)
# scatter!(0, 0, alpha = 0.1, color = :red)
