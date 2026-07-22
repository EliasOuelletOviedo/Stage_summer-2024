using FileIO
using ImageIO
using CairoMakie
using ProgressBars

# folder = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Test galvo/60 Hz x 600 images-001/Bliq VMS/C2"
folder = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Test galvo/1024x512 images-007/Bliq VMS/C2"

tiff_files = filter(f ->
    endswith(lowercase(f), ".tif") ||
    endswith(lowercase(f), ".tiff"),
    readdir(folder, join=true)
)

# Fonction qui extrait le numéro dans "T123.tif"
function frame_number(filepath)
    filename = basename(filepath)
    m = match(r"T(\d+)\.tif{1,2}f?$"i, filename)
    return parse(Int, m.captures[1])
end

sort!(tiff_files, by=frame_number)

first_profile = vec(sum(Float64.(load(tiff_files[1])), dims=1))

n_images = length(tiff_files)
profile_length = length(first_profile)

# Chaque ligne = un fichier TIFF
big_image = Matrix{Float64}(undef, n_images, profile_length)

for (i, file) in ProgressBar(enumerate(tiff_files))
    data = Float64.(load(file))
    big_image[i, :] .= vec(sum(data, dims=2))
end

println(size(big_image))

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect(), yreversed = true)

image!(ax, log10.(big_image[86:506,:]' .+ 1))
# image!(ax, log10.(big_image' .+ 1))
# lines!(ax, big_image[560, :], color=:red)

ax.aspect = 1

fig