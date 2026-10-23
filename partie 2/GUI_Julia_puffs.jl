using LibSerialPort
using GLMakie
using Observables

set_theme!(backgroundcolor = :gray95)

# Set colors
dark_gray = RGBf(0.4, 0.4, 0.4)
mid_gray = RGBf(0.6, 0.6, 0.6)
light_gray = RGBf(0.8, 0.8, 0.8)
off_white = RGBf(0.9, 0.9, 0.9)

fig = Figure(size = (1000, 600))

fig[1, 1] = plotGrid = GridLayout(tellwidth = true)
fig[2, 1] = controlGrid = GridLayout(tellwidth = false)
controlGrid[1, 1] = buttonGrid = GridLayout(tellwidth = false)
controlGrid[1, 2] = sliderGrid = GridLayout(tellwidth = false)
buttonGrid[2, 1] = pulsewidthGrid = GridLayout(tellwidth = false)
buttonGrid[2, 2] = frequencyGrid  = GridLayout(tellwidth = false)

global input_axis = plotGrid[1, 1] = Axis(fig,
    title = "Input",
    ylabel = "Intensity [%]",
    titlesize = 16,
    aspect = 3,
    xtickformat = "{:d}")

global sigmoid_axis = plotGrid[2, 1] = Axis(fig,
    title = "Sigmoid",
    ylabel = "Intensity [%]",
    titlesize = 16,
    aspect = 3,
    limits = (-5.25, 5.25, -0.05, 1.05))

global analog_axis = plotGrid[1, 2] = Axis(fig,
    title = "Analog output",
    ylabel = "Intensity [%]",
    titlesize = 16,
    aspect = 3,
    xtickformat = "{:d}")

global digital_axis = plotGrid[2, 2] = Axis(fig,
    title = "Digital output",
    ylabel = "Intensity [%]",
    titlesize = 16,
    aspect = 3,
    xtickformat = "{:d}")

pulsewidth_label = Label(buttonGrid[1, 1],
fontsize = 16,
justification = :right,
text = "Pusle width [ms] :")

pulsewidth_textbox = Textbox(pulsewidthGrid[1:2, 1],
    borderwidth = 0,
    boxcolor = :white,
    boxcolor_focused = :white,
    boxcolor_focused_invalid = RGBf(1, 0.9, 0.9),
    boxcolor_hover = off_white,
    cornerradius = 0,
    height = 32,
    placeholder = " ",
    validator = Float64,
    width = 96)
pulsewidth_textbox.displayed_string[] = "10.0"
pulsewidth_textbox.stored_string[] = "10.0"

# Frequency Up/Down Buttons
pulsewidth_up = Button(pulsewidthGrid[1, 2],
    buttoncolor = :gray95,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = light_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 8,
    height = 17,
    label = "▲",
    labelcolor_active = :black,
    strokecolor = dark_gray,
    strokewidth = 1,
    valign = :bottom,
    width = 32)

pulsewidth_down = Button(pulsewidthGrid[2, 2],
    buttoncolor = :gray95,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = light_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 8,
    height = 17,
    label = "▼",
    labelcolor_active = :black,
    strokecolor = dark_gray,
    strokewidth = 1,
    valign = :top,
    width = 32)

# A transparent box for bordering
Box(pulsewidthGrid[1:2, 1:2],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)

frequency_label = Label(buttonGrid[1, 2],
fontsize = 16,
justification = :right,
text = "Image frequency [ms] :")

frequency_textbox = Textbox(frequencyGrid[1:2, 1],
    borderwidth = 0,
    boxcolor = :white,
    boxcolor_focused = :white,
    boxcolor_focused_invalid = RGBf(1, 0.9, 0.9),
    boxcolor_hover = off_white,
    cornerradius = 0,
    height = 32,
    placeholder = " ",
    validator = Float64,
    width = 96)
frequency_textbox.displayed_string[] = "10.0"
frequency_textbox.stored_string[] = "10.0"

# Frequency Up/Down Buttons
frequency_up = Button(frequencyGrid[1, 2],
    buttoncolor = :gray95,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = light_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 8,
    height = 17,
    label = "▲",
    labelcolor_active = :black,
    strokecolor = dark_gray,
    strokewidth = 1,
    valign = :bottom,
    width = 32)

frequency_down = Button(frequencyGrid[2, 2],
    buttoncolor = :gray95,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = light_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 8,
    height = 17,
    label = "▼",
    labelcolor_active = :black,
    strokecolor = dark_gray,
    strokewidth = 1,
    valign = :top,
    width = 32)

# A transparent box for bordering
Box(frequencyGrid[1:2, 1:2],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)


# Set gap adjustments
colgap!(pulsewidthGrid, 1, Fixed(-1))
rowgap!(pulsewidthGrid, 1, Fixed(0))
colgap!(frequencyGrid, 1, Fixed(-1))
rowgap!(frequencyGrid, 1, Fixed(0))
# colgap!(fpsGrid, 1, Fixed(-1))
# rowgap!(fpsGrid, 1, Fixed(0))
colgap!(buttonGrid, 1, Fixed(64))

params_sliders = SliderGrid(
    sliderGrid[1, 1],
    (label = "Proportional", range = 0:1:100, format = "{:d}", startvalue = 0, linewidth = 14,
        color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
    (label = "Integral", range = 0:1:100, format = "{:d}", startvalue = 0, linewidth = 14,
        color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
    (label = "Threshold", range = 0:1:100, format = "{:d}", startvalue = 50, linewidth = 14,
        color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
    width = 350)

# params_sliders = SliderGrid(
#     sliderGrid[1, 1],
#     (label = "Sigmoid width", range = -5:0.1:5, format = "{:.1f}", startvalue = 0, linewidth = 14,
#         color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
#     (label = "Sigmoid offset", range = -5:0.1:5, format = "{:.1f}", startvalue = 0, linewidth = 14,
#         color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
#     (label = "Threshold", range = 0:1:100, format = "{:d}", startvalue = 50, linewidth = 14,
#         color_active=dark_gray, color_active_dimmed=mid_gray, color_inactive=light_gray, tellheight=false),
#     width = 350)

t = Observable([0])
input = Observable([0])
ana_out = Observable([0])
dig_out = Observable([0])

base_x = range(-5, 5, 1001)

x = Observable(base_x .* 10^(params_sliders.sliders[1].value[]))

function sig_func(x, a, b)
    return 1 / (1 + exp(-10^a * (x + b)))
end

sigmoid = sig_func.(base_x, 0, 0)

lines!(sigmoid_axis, x, sigmoid, color=Makie.wong_colors()[1])

on(params_sliders.sliders[1].value) do val
    x[] = (base_x .- params_sliders.sliders[2].value[]) .* 10^val
    sigmoid = sig_func.(base_x, val, params_sliders.sliders[2].value[])

    sigmoid_axis.limits[] = ((-5.25 - params_sliders.sliders[2].value[]) * 10^val, (5.25 - params_sliders.sliders[2].value[]) * 10^val, -0.05, 1.05)
end

on(params_sliders.sliders[2].value) do val
    x[] = (base_x .- val) .* 10^params_sliders.sliders[1].value[]
    sigmoid = sig_func.(base_x, params_sliders.sliders[1].value[], val)

    sigmoid_axis.limits[] = ((-5.25 - val) * 10^params_sliders.sliders[1].value[], (5.25 - val) * 10^params_sliders.sliders[1].value[], -0.05, 1.05)
end

fig
