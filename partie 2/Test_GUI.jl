using GLMakie

# 1. List of common system font families
fonts = [
    "Arial",
    "Helvetica",      # often aliased to Arial on some platforms
    "DejaVu Sans",
    "Liberation Sans",
    "Roboto",
    "Open Sans",
    "Times New Roman",
    "Georgia",
    "DejaVu Serif",
    "Liberation Serif",
    "Noto Serif",
    "Playfair Display",
    "Courier New",
    "Consolas",
    "DejaVu Sans Mono",
    "Liberation Mono",
    "Source Code Pro",
    "Fira Code",
]

set_theme!(backgroundcolor = :gray12, textcolor = :gray85)

# 2. Create a tall, scrollable figure
n = length(fonts)
fig = Figure(size = (1440, 847))

println(n)

global row = 1

# 3. For each font, add a label
for (i, fam) in enumerate(fonts)
    ax = Axis(fig[row, i + 3 - 3*row];
        xgridvisible = false, ygridvisible = false,
        title = "", backgroundcolor = :gray8
    )
    hidedecorations!(ax)

    # Write the name in its own font
    text!(ax, 0.5, 0.5,
        text = fam,
        align = (:center, :center),
        font = fam,
        fontsize = 24,
    )
    hidespines!(ax)

    if i == 3
        global row = row + 1
    elseif i == 6
        global row = row + 1
    elseif i == 12
        global row = row + 1
    elseif i == 15
        global row = row + 1
    elseif i == 9
        global row = row + 1
    end
end

# 4. Display
fig
