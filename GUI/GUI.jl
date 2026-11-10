using GLMakie
using FFTW
using Statistics
using NativeFileDialog
using ZipFile
using Observables
using Dates
using LibSerialPort
using CSV
using DataFrames
using GLMakie.GLFW
using GLMakie: to_native

include("functions.jl")
include("pop_up.jl")
include("pop_up_2.jl")

portname = "usbmodem102"

try
    LibSerialPort.close(serial)
catch
end

# Common plot and widget style attributes
function range_validator(a::Integer, b::Integer)
    """
    Returns a validator function ensuring parsed integer s is in [a, b].
    """
    @assert a ≤ b

    return function(s::String)
        s_trim = strip(s)
        x = tryparse(Int, s_trim)
        return x !== nothing && a ≤ x ≤ b
    end
end

BUTTON_ATTRS       = Dict(
    :buttoncolor         => light_gray,
    :buttoncolor_active  => dark_gray,
    :buttoncolor_hover   => mid_gray,
    :cornerradius        => 4,
    :fontsize            => 16,
    :height              => 32,
    :labelcolor_hover    => :white,
    :width               => 112
)

MENU_ATTRS         = Dict(
    :cell_color_active             => mid_gray,
    :cell_color_inactive_even      => light_gray,
    :cell_color_inactive_odd       => light_gray,
    :cell_color_hover              => mid_gray,
    :fontsize                      => 16,
    :height                        => 32,
    :selection_cell_color_inactive => light_gray,
    :textcolor                     => :gray90,
    :width                         => 112
)

SPINNER_BOX_ATTRS  = Dict(
    :color        => mid_gray,
    :cornerradius => 2,
    :height       => 24,
    :width        => 88,
    :strokecolor  => mid_gray,
    :strokewidth  => 0.5,
)

SPINNER_TEXT_ATTRS = Dict(
    :borderwidth                 => 0,
    :boxcolor                    => dark_gray,
    :boxcolor_focused            => mid_gray,
    :boxcolor_focused_invalid    => RGBf(0.22, 0.11, 0.11),
    :boxcolor_hover              => light_gray,
    :cornerradius                => 0,
    :halign                      => :left,
    :height                      => 24,
    :placeholder                 => " ",
    :reset_on_defocus            => true,
    :textcolor_placeholder       => :black,
    :textpadding                 => (8, 0, 0, 4),
    :width                       => 64
)

SPINNER_UP_ATTRS   = Dict(
    :buttoncolor        => mid_gray,
    :buttoncolor_active => dark_gray,
    :buttoncolor_hover  => light_gray,
    :cornerradius       => 0,
    :font               => :bold,
    :fontsize           => 6,
    :halign             => :right,
    :height             => 12,
    :label              => "▲",
    :labelcolor_active  => :white,
    :labelcolor_hover   => :white,
    :strokewidth        => 0,
    :valign             => :top,
    :width              => 24
)

SPINNER_DOWN_ATTRS = Dict(
    :buttoncolor        => mid_gray,
    :buttoncolor_active => dark_gray,
    :buttoncolor_hover  => light_gray,
    :cornerradius       => 0,
    :font               => :bold,
    :fontsize           => 6,
    :halign             => :right,
    :height             => 12,
    :label              => "▼",
    :labelcolor_active  => :white,
    :labelcolor_hover   => :white,
    :strokewidth        => 0,
    :valign             => :bottom,
    :width              => 24
)

# Build main figure and layout
fig = Figure(size = (1440, 847), msaa = 1)
plotgrid     = GridLayout(fig[1, 1],          tellwidth=true)
buttongrid   = GridLayout(fig[2, 1],          tellwidth=false)
progressgrid = GridLayout(fig[3, 1],          tellwidth=false)
textboxgrid  = GridLayout(buttongrid[1:3, 4], tellwidth=false)
paramsgrid   = GridLayout(buttongrid[1:3, 7], tellwidth=true)

# Create Axes
global axis_hist   = Axis(plotgrid[1, 1]; merge(AXIS_ATTRS, Dict(
    :title  => "Histogram",
    :xlabel => "Time [ns]",
    :ylabel => "Photon counts"
))...)

global axis_counts = Axis(plotgrid[2, 1]; merge(AXIS_ATTRS, Dict(
    :title  => "Photon counts",
    :xlabel => "Time [s]",
    :ylabel => "Photon counts"
))...)

global axis_life   = Axis(plotgrid[1, 2]; merge(AXIS_ATTRS, Dict(
    # :title  => "Ion concentration",
    :title  => "Lifetime",
    :xlabel => "Time [s]", 
    # :ylabel => "Ion concentration [mM]"
    :ylabel => "Lifetime [ns]"
))...)

global axis_signal = Axis(plotgrid[2, 2]; merge(AXIS_ATTRS, Dict(
    :title  => "Light intensity",
    :xlabel => "Time [s]",
    :ylabel => "Duty cycle [%]"
))...)

ylims!(axis_signal, -5, 105)

# Labels and Controls
label_frequency   = Label(buttongrid[1, 1],       text="Frequency : 0 Hz", width = 120)
label_number      = Label(buttongrid[2, 1],       text="Files : 0",        width = nothing)
label_inverted    = Label(buttongrid[3, 1][1, 1], text="Inverted",         width = nothing)
checkbox_inverted = Toggle(
    buttongrid[3, 1][1, 2],
    active              = true,
    rimfraction         = 0.25,
    buttoncolor         = dark_gray,
    framecolor_active   = mid_gray,
    framecolor_inactive = light_gray
)

# Buttons
button_start    = Button(buttongrid[1, 2]; merge(BUTTON_ATTRS, Dict(:label => "START"))...)
button_stop     = Button(buttongrid[1, 3]; merge(BUTTON_ATTRS, Dict(:label => "CLEAR"))...)
button_load     = Button(buttongrid[2, 2]; merge(BUTTON_ATTRS, Dict(:label => "LOAD"))...)
button_load_IRF = Button(buttongrid[2, 3]; merge(BUTTON_ATTRS, Dict(:label => "LOAD IRF"))...)
button_pop_up   = Button(buttongrid[3, 2]; merge(BUTTON_ATTRS, Dict(:label => "PROTOCOL"))...)
menu_mode       = Menu(buttongrid[3, 3];   merge(MENU_ATTRS,   Dict(
    :default => "Playback",
    :options => ["Analysis", "Playback", "Real-time", "Save"]
))...)

# Textboxes & Labels for serial & timing parameters
portname_textbox = Textbox(textboxgrid[1, 1:2]; merge(TEXTBOX_ATTRS, Dict(
    :displayed_string => portname,
    :stored_string    => portname,
    :width            => nothing
))...)
textbox_fps      = Textbox(textboxgrid[2, 2];   merge(TEXTBOX_ATTRS, Dict(
    :displayed_string => "50",
    :stored_string    => "50"
))...)
textbox_refresh  = Textbox(textboxgrid[3, 2];   merge(TEXTBOX_ATTRS, Dict(
    :displayed_string => "5",
    :stored_string    => "5"
))...)
textbox_min      = Textbox(textboxgrid[4, 2];   merge(TEXTBOX_ATTRS, Dict(
    :displayed_string => "1",
    :stored_string    => "1"
))...)
label_fps        = Label(textboxgrid[2, 1], text = "Blue light [%] :", halign = :right, width = nothing)
label_min        = Label(textboxgrid[3, 1], text = "Refresh [Hz] :",   halign = :right, width = nothing)
label_max        = Label(textboxgrid[4, 1], text = "Maximum [%] :",    halign = :right, width = nothing)

# Activate/connect controls
button_active  = Button(buttongrid[1, 5]; merge(BUTTON_ATTRS, Dict(
    :label             => "Activate",
    :buttoncolor       => mid_gray,
    :buttoncolor_hover => dark_gray
))...)
button_connect = Button(buttongrid[2, 5]; merge(BUTTON_ATTRS, Dict(
    :label             => "Connect",
    :buttoncolor       => mid_gray,
    :buttoncolor_hover => dark_gray
))...)
menu_output    = Menu(buttongrid[3, 5];   merge(MENU_ATTRS,   Dict(
    :default => "Digital",
    :options => ["Analog","Digital"]
))...)

# Monitor box
monitor_box = Box(buttongrid[1:3, 6],
    color = dark_gray,
    cornerradius = 2,
    height = nothing,
    strokewidth = 0.5,
    strokecolor = light_gray,
    width = 210
)
global monitor = Label(buttongrid[1:3, 6],
    "\n",
    fontsize = 15,
    height = nothing,
    justification = :left,
    padding = (8, 8, 8, 8),
    tellheight = false,
    tellwidth = false,
    valign = :bottom,
    width = 210
)

# Upper parameters grid
Label(paramsgrid[1, 1:4]; text="General settings", halign=:center, height=nothing, tellwidth=true)

Label(paramsgrid[2, 1]; text="Time range", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[3, 1];   SPINNER_BOX_ATTRS...)
range_up      = Button(paramsgrid[3, 1];  SPINNER_UP_ATTRS...)
range_down    = Button(paramsgrid[3, 1];  SPINNER_DOWN_ATTRS...)
range_textbox = Textbox(paramsgrid[3, 1]; merge(SPINNER_TEXT_ATTRS, Dict(
    :displayed_string => "600",
    :stored_string    => "600",
    :validator        => UInt32
))...)

Label(paramsgrid[2, 2]; text="Binning", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[3, 2];   SPINNER_BOX_ATTRS...)
averaging_up      = Button(paramsgrid[3, 2];  SPINNER_UP_ATTRS...)
averaging_down    = Button(paramsgrid[3, 2];  SPINNER_DOWN_ATTRS...)
averaging_textbox = Textbox(paramsgrid[3, 2]; merge(SPINNER_TEXT_ATTRS, Dict(
    :displayed_string => "1",
    :stored_string    => "1",
    :validator        => range_validator(1,100)
))...)

Label(paramsgrid[2, 3]; text="Lifetimes", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[3, 3];   SPINNER_BOX_ATTRS...)
lifetimes_up      = Button(paramsgrid[3, 3];  SPINNER_UP_ATTRS...)
lifetimes_down    = Button(paramsgrid[3, 3];  SPINNER_DOWN_ATTRS...)
lifetimes_textbox = Textbox(paramsgrid[3, 3]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "1",
    :stored_string    => "1",
    :validator        => range_validator(1,3)
))...)

Label(paramsgrid[2, 4]; text="Smoothing", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[3, 4];   SPINNER_BOX_ATTRS...)
smoothing_up      = Button(paramsgrid[3, 4];  SPINNER_UP_ATTRS...)
smoothing_down    = Button(paramsgrid[3, 4];  SPINNER_DOWN_ATTRS...)
smoothing_textbox = Textbox(paramsgrid[3, 4]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "5",
    :stored_string    => "5",
    :validator        => range_validator(0,10)
))...)

# Lower parameters grid
Label(paramsgrid[5, 1:4]; text="Controller settings", halign=:center, height=nothing, tellwidth=true)

Label(paramsgrid[6, 1]; text="Proportionnal", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[7, 1];   SPINNER_BOX_ATTRS...)
proportionnal_up      = Button(paramsgrid[7, 1];  SPINNER_UP_ATTRS...)
proportionnal_down    = Button(paramsgrid[7, 1];  SPINNER_DOWN_ATTRS...)
proportionnal_textbox = Textbox(paramsgrid[7, 1]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "0",
    :stored_string    => "0",
    :validator        => UInt32
))...)

Label(paramsgrid[6, 2]; text="Integral", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[7, 2];   SPINNER_BOX_ATTRS...)
integral_up      = Button(paramsgrid[7, 2];  SPINNER_UP_ATTRS...)
integral_down    = Button(paramsgrid[7, 2];  SPINNER_DOWN_ATTRS...)
integral_textbox = Textbox(paramsgrid[7, 2]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "0",
    :stored_string    => "0",
    :validator        => UInt32
))...)

Label(paramsgrid[6, 3]; text="Derivative", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[7, 3];   SPINNER_BOX_ATTRS...)
derivative_up      = Button(paramsgrid[7, 3];  SPINNER_UP_ATTRS...)
derivative_down    = Button(paramsgrid[7,  3];  SPINNER_DOWN_ATTRS...)
derivative_textbox = Textbox(paramsgrid[7  ,3]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "0",
    :stored_string    => "0",
    :validator        => Float64
))...) 

Label(paramsgrid[6, 4]; text="Setpoint", halign=:center, height=nothing, tellwidth=true)
Box(paramsgrid[7, 4];   SPINNER_BOX_ATTRS...)
setpoint_up      = Button(paramsgrid[7, 4];  SPINNER_UP_ATTRS...)
setpoint_down    = Button(paramsgrid[7, 4];  SPINNER_DOWN_ATTRS...)
setpoint_textbox = Textbox(paramsgrid[7, 4]; merge(SPINNER_TEXT_ATTRS,Dict(
    :displayed_string => "5.0",
    :stored_string    => "5.0",
    :validator        => Float64
))...)

# Progress bar
Box(progressgrid[1, 1],
    color        = light_gray,
    cornerradius = 5,
    halign       = :left,
    height       = 10,
    strokecolor  = light_gray,
    strokewidth  = 1,
    width        = 1300
)
global box_progress   = Box(progressgrid[1, 1],
    color         = dark_gray,
    cornerradius  = 5,
    halign        = :left,
    height        = 10,
    strokevisible = false,
    tellheight    = false,
    tellwidth     = false,
    width         = 0
)
global label_progress = Label(progressgrid[1, 3], text = "0 %", halign = :right)

# Layout adjustments
colsize!(plotgrid,     1, 640)
colsize!(plotgrid,     2, 640)
colsize!(buttongrid,   4, 160)
colsize!(paramsgrid,   1, 88)
colsize!(paramsgrid,   2, 88)
colsize!(paramsgrid,   3, 88)
colsize!(paramsgrid,   4, 88)
colsize!(progressgrid, 2, 10)
colsize!(progressgrid, 3, 40)
colgap!(buttongrid,    3, 30)
# rowsize!(paramsgrid,   1, 10)
colgap!(fig.layout,    0)
colgap!(buttongrid,    20)
colgap!(textboxgrid,   10)
colgap!(progressgrid,  0)
rowgap!(buttongrid,    20)
rowgap!(fig.layout,    1, 40)
rowgap!(textboxgrid,   12)
rowgap!(paramsgrid,    2, 4)
rowgap!(paramsgrid,    6, 4)

# Global data & state
global irf               = get_irf()
global irf_bin_size      = get_irf_bin_size()
global tcspc_window_size = round(irf[end, 1] + irf[2, 1], sigdigits=4)
global fft_plan          = plan_fft(zeros(Float64, 256))
global ifft_plan         = plan_ifft(zeros(Float64, 256))
global folderpath        = get_folder_path("FLIM")
global serial            = nothing
global running           = Observable(false)
global pause             = Observable(false)
global active            = Observable(false)
global connected         = Observable(false)
global reset_I           = Observable(false)
global protocol_active   = Observable(false)

setpoint = hlines!(axis_life, parse(Float64, setpoint_textbox.stored_string[]), color=Makie.wong_colors()[6])

# Event handlers
on(button_start.clicks) do _
    if !running[]
        running[] = true

        button_start.label = "PAUSE"
        button_stop.label = "STOP"

        box_progress.width = 0
        label_progress.text = "0 %"
        label_frequency.text = "Frequency : 0 Hz"
        label_number.text = "Files : 0"

        empty!(axis_hist)
        empty!(axis_counts)
        empty!(axis_life)
        empty!(axis_signal)

        global setpoint = hlines!(axis_life, parse(Float64, setpoint_textbox.stored_string[]), color=Makie.wong_colors()[6])

        if menu_mode.selection[] == "Playback"
            if protocol_active[]
                @async start_protocol_playback()
            else
                @async start_playback()
            end
            # @async start_protocol_playback()
        elseif menu_mode.selection[] == "Analysis"
            start_analysis()
        elseif menu_mode.selection[] == "Save"
            @async start_save()
        elseif menu_mode.selection[] == "Real-time"
            if protocol_active[]
                @async start_protocol_real_time()
            else
                @async start_real_time()
            end
        end
    elseif pause[]
        pause[] = false

        button_start.label = "PAUSE"
        button_stop.label = "STOP"
    else
        pause[] = true

        button_start.label = "START"
        button_stop.label = "CLEAR"
    end
end

on(button_stop.clicks) do _
    running[] = false
    pause[] = false

    if connected[]
        # clear(serial)
        send_command(serial, "A 0 DS 1 0\n")
    end

    button_start.label = "START"
    button_stop.label = "CLEAR"

    box_progress.width = 0
    label_progress.text = "0 %"
    label_frequency.text = "Frequency : 0 Hz"
    label_number.text = "Files : 0"

    empty!(axis_hist)
    empty!(axis_counts)
    empty!(axis_life)
    empty!(axis_signal)

    global setpoint = hlines!(axis_life, parse(Float64, setpoint_textbox.stored_string[]), color=Makie.wong_colors()[6])
end

on(button_load.clicks) do _
    global folderpath = get_new_folder_path("FLIM")
end

on(button_load_IRF.clicks) do _
    global irf = get_new_irf()
    global irf_bin_size = get_irf_bin_size()
end

on(button_pop_up.clicks) do _
    try
        if pop_up_screen.window_open[]
            GLFW.SetWindowShouldClose(pop_up_window, true)
        end
        sleep(0.1)
        global pop_up_screen = display(GLMakie.Screen(), fig_2)
        global pop_up_window = to_native(pop_up_screen)
    catch e
        global pop_up_screen = display(GLMakie.Screen(), fig_2)
        global pop_up_window = to_native(pop_up_screen)
    end
    sleep(0.1)
    try
        if pop_up_2_screen.window_open[]
            GLFW.SetWindowShouldClose(pop_up_2_window, true)
        end
        sleep(0.1)
        global pop_up_2_screen = display(GLMakie.Screen(), fig_3)
        global pop_up_2_window = to_native(pop_up_2_screen)
    catch e
        global pop_up_2_screen = display(GLMakie.Screen(), fig_3)
        global pop_up_2_window = to_native(pop_up_2_screen)
    end
    sleep(1)
end

on(button_connect.clicks) do _
    if connected[]
        clear(serial)
        LibSerialPort.close(serial)
        add_to_monitor("Disconnection successful")
        add_to_monitor("Deactivated")
        println("Disconnection successful")

        button_connect.label[] = "Connect"
        button_active.label = "Activate"
        connected[] = false
        active[] = false

        button_active.buttoncolor[] = mid_gray
        button_active.buttoncolor_hover[] = dark_gray

        button_connect.buttoncolor[] = mid_gray
        button_connect.buttoncolor_hover[] = dark_gray
    else
        global serial = connect_to_port(portname_textbox.stored_string[])

        sleep(0.1)
        
        try
            send_command(serial, "C Z\nC R\n")
            send_command(serial, "S 1 S 1\n")

            button_connect.label[] = "Disconnect"

            button_connect.buttoncolor[] = light_gray
            button_connect.buttoncolor_hover[] = mid_gray
        catch e
            println(e)
        end
    end
end

on(button_active.clicks) do _
    if connected[]
        if active[]
            active[] = false

            button_active.label = "Activate"
            button_active.buttoncolor[] = mid_gray
            button_active.buttoncolor_hover[] = dark_gray

            add_to_monitor("Deactivated")
            clear(serial)
        else
            active[] = true

            button_active.label = "Deactivate"
            button_active.buttoncolor[] = light_gray
            button_active.buttoncolor_hover[] = mid_gray

            add_to_monitor("Activated")
            
            @async start_signal(
                parse(Float64, textbox_fps.stored_string[]),
                parse(Float64, textbox_refresh.stored_string[]),
                parse(Float64, textbox_min.stored_string[])
            )
        end
    end
end

on(range_up.clicks) do _
    if !isnothing(tryparse(UInt32, range_textbox.displayed_string[]))
        value = parse(UInt32, range_textbox.displayed_string[])
        multiplier = floor(Int, round(log(10, value), digits = 4))
        res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

        range_textbox.stored_string[] = string(res)
        range_textbox.displayed_string[] = string(res)
    else
        range_textbox.stored_string[] = "60"
        range_textbox.displayed_string[] = "60"
    end
end

on(range_down.clicks) do _
    if !isnothing(tryparse(UInt32, range_textbox.displayed_string[]))
        value = parse(UInt32, range_textbox.displayed_string[])

        if value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            range_textbox.stored_string[] = string(res)
            range_textbox.displayed_string[] = string(res)
        end
    else
        range_textbox.stored_string[] = "60"
        range_textbox.displayed_string[] = "60"
    end
end

on(averaging_up.clicks) do _
    if !isnothing(tryparse(UInt32, averaging_textbox.displayed_string[]))
        value = parse(UInt32, averaging_textbox.displayed_string[])

        if value < 1 | value > 100
            averaging_textbox.stored_string[] = "1"
            averaging_textbox.displayed_string[] = "1"
        end
        
        if value != 100
            averaging_textbox.stored_string[] = string(value + 1)
            averaging_textbox.displayed_string[] = string(value + 1)
        end
    else
        averaging_textbox.stored_string[] = "1"
        averaging_textbox.displayed_string[] = "1"
    end
end

on(averaging_down.clicks) do _
    if !isnothing(tryparse(UInt32, averaging_textbox.displayed_string[]))
        value = parse(UInt32, averaging_textbox.displayed_string[])

        if value < 1 | value > 100
            averaging_textbox.stored_string[] = "1"
            averaging_textbox.displayed_string[] = "1"
        end

        if value != 1
            averaging_textbox.stored_string[] = string(value - 1)
            averaging_textbox.displayed_string[] = string(value - 1)
        end
    else
        averaging_textbox.stored_string[] = "1"
        averaging_textbox.displayed_string[] = "1"
    end
end

on(lifetimes_up.clicks) do _
    if !isnothing(tryparse(UInt32, lifetimes_textbox.displayed_string[]))
        value = parse(UInt32, lifetimes_textbox.displayed_string[])

        if value == 1
            lifetimes_textbox.stored_string[] = "2"
            lifetimes_textbox.displayed_string[] = "2"
        elseif value == 2
            lifetimes_textbox.stored_string[] = "3"
            lifetimes_textbox.displayed_string[] = "3"
        elseif value == 3
        else 
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        end
    else
        lifetimes_textbox.stored_string[] = "1"
        lifetimes_textbox.displayed_string[] = "1"
    end
end

on(lifetimes_down.clicks) do _
    if !isnothing(tryparse(UInt32, lifetimes_textbox.displayed_string[]))
        value = parse(UInt32, lifetimes_textbox.displayed_string[])

        if value == 1
        elseif value == 2
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        elseif value == 3
            lifetimes_textbox.stored_string[] = "2"
            lifetimes_textbox.displayed_string[] = "2"
        else 
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        end
    else
        lifetimes_textbox.stored_string[] = "1"
        lifetimes_textbox.displayed_string[] = "1"
    end
end

on(smoothing_up.clicks) do _
    if !isnothing(tryparse(UInt32, smoothing_textbox.displayed_string[]))
        value = parse(UInt32, smoothing_textbox.displayed_string[])

        if value < 0 | value > 5
            smoothing_textbox.stored_string[] = "0"
            smoothing_textbox.displayed_string[] = "0"
        end
        # println(value)
        if value != 10
            smoothing_textbox.stored_string[] = string(value + 1)
            smoothing_textbox.displayed_string[] = string(value + 1)
        end
    else
        smoothing_textbox.stored_string[] = "0"
        smoothing_textbox.displayed_string[] = "0"
    end
end

on(smoothing_down.clicks) do _
    if !isnothing(tryparse(UInt32, smoothing_textbox.displayed_string[]))
        value = parse(UInt32, smoothing_textbox.displayed_string[])

        if value < 0 | value > 10
            smoothing_textbox.stored_string[] = "0"
            smoothing_textbox.displayed_string[] = "0"
        end

        if value != 0
            smoothing_textbox.stored_string[] = string(value - 1)
            smoothing_textbox.displayed_string[] = string(value - 1)
        end
    else
        smoothing_textbox.stored_string[] = "0"
        smoothing_textbox.displayed_string[] = "0"
    end
end

# on(smoothing_up.clicks) do _
#     if !isnothing(tryparse(UInt32, smoothing_textbox.displayed_string[]))
#         value = parse(UInt32, smoothing_textbox.displayed_string[])
#         multiplier = floor(Int, round(log(10, value), digits = 4))
#         res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

#         smoothing_textbox.stored_string[] = string(res)
#         smoothing_textbox.displayed_string[] = string(res)
#     else
#         smoothing_textbox.stored_string[] = "60"
#         smoothing_textbox.displayed_string[] = "60"
#     end
# end

# on(smoothing_down.clicks) do _
#     if !isnothing(tryparse(UInt32, smoothing_textbox.displayed_string[]))
#         value = parse(UInt32, smoothing_textbox.displayed_string[])

#         if value != 1
#             multiplier = floor(Int, round(log(10, value - 1), digits = 8))
#             res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

#             smoothing_textbox.stored_string[] = string(res)
#             smoothing_textbox.displayed_string[] = string(res)
#         end
#     else
#         smoothing_textbox.stored_string[] = "60"
#         smoothing_textbox.displayed_string[] = "60"
#     end
# end

on(proportionnal_up.clicks) do _
    if !isnothing(tryparse(UInt32, proportionnal_textbox.displayed_string[]))
        value = parse(UInt32, proportionnal_textbox.displayed_string[])
        if value != 0
            multiplier = floor(Int, round(log(10, value), digits = 4))
            res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

            proportionnal_textbox.stored_string[] = string(res)
            proportionnal_textbox.displayed_string[] = string(res)
        else
            proportionnal_textbox.stored_string[] = "1"
            proportionnal_textbox.displayed_string[] = "1"
        end
    else
        proportionnal_textbox.stored_string[] = "0"
        proportionnal_textbox.displayed_string[] = "0"
    end
end

on(proportionnal_down.clicks) do _
    if !isnothing(tryparse(UInt32, proportionnal_textbox.displayed_string[]))
        value = parse(UInt32, proportionnal_textbox.displayed_string[])

        if value == 0
        elseif value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            proportionnal_textbox.stored_string[] = string(res)
            proportionnal_textbox.displayed_string[] = string(res)
        else
            proportionnal_textbox.stored_string[] = "0"
            proportionnal_textbox.displayed_string[] = "0"
        end
    else
        proportionnal_textbox.stored_string[] = "0"
        proportionnal_textbox.displayed_string[] = "0"
    end
end

on(integral_up.clicks) do _
    if !isnothing(tryparse(UInt32, integral_textbox.displayed_string[]))
        value = parse(UInt32, integral_textbox.displayed_string[])
        if value != 0
            multiplier = floor(Int, round(log(10, value), digits = 4))
            res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

            integral_textbox.stored_string[] = string(res)
            integral_textbox.displayed_string[] = string(res)
        else
            integral_textbox.stored_string[] = "1"
            integral_textbox.displayed_string[] = "1"
        end
    else
        integral_textbox.stored_string[] = "0"
        integral_textbox.displayed_string[] = "0"
    end
end

on(integral_down.clicks) do _
    if !isnothing(tryparse(UInt32, integral_textbox.displayed_string[]))
        value = parse(UInt32, integral_textbox.displayed_string[])

        if value == 0
        elseif value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            integral_textbox.stored_string[] = string(res)
            integral_textbox.displayed_string[] = string(res)
        else
            integral_textbox.stored_string[] = "0"
            integral_textbox.displayed_string[] = "0"
        end
    else
        integral_textbox.stored_string[] = "0"
        integral_textbox.displayed_string[] = "0"
    end
end

on(derivative_up.clicks) do _
    if !isnothing(tryparse(UInt32, derivative_textbox.displayed_string[]))
        value = parse(UInt32, derivative_textbox.displayed_string[])
        if value != 0
            multiplier = floor(Int, round(log(10, value), digits = 4))
            res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

            derivative_textbox.stored_string[] = string(res)
            derivative_textbox.displayed_string[] = string(res)
        else
            derivative_textbox.stored_string[] = "1"
            derivative_textbox.displayed_string[] = "1"
        end
    else
        derivative_textbox.stored_string[] = "0"
        derivative_textbox.displayed_string[] = "0"
    end
end

on(derivative_down.clicks) do _
    if !isnothing(tryparse(UInt32, derivative_textbox.displayed_string[]))
        value = parse(UInt32, derivative_textbox.displayed_string[])

        if value == 0
        elseif value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            derivative_textbox.stored_string[] = string(res)
            derivative_textbox.displayed_string[] = string(res)
        else
            derivative_textbox.stored_string[] = "0"
            derivative_textbox.displayed_string[] = "0"
        end
    else
        derivative_textbox.stored_string[] = "0"
        derivative_textbox.displayed_string[] = "0"
    end
end

on(setpoint_up.clicks) do _
    if !isnothing(tryparse(Float64, setpoint_textbox.displayed_string[]))
        value = parse(Float64, setpoint_textbox.displayed_string[])
        if value < 0
            setpoint_textbox.stored_string[] = "0"
            setpoint_textbox.displayed_string[] = "0"
        elseif value == 0
            setpoint_textbox.stored_string[] = "0.1"
            setpoint_textbox.displayed_string[] = "0.1"
        else
            multiplier = max(floor(Int, round(log(10, value) - 1, digits = 4)), -1)

            setpoint_textbox.stored_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
            setpoint_textbox.displayed_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
        end
    else
        setpoint_textbox.stored_string[] = "5"
        setpoint_textbox.displayed_string[] = "5"
    end
end

on(setpoint_down.clicks) do _
    if !isnothing(tryparse(Float64, setpoint_textbox.displayed_string[]))
        value = parse(Float64, setpoint_textbox.displayed_string[])
        if value <= 0
            setpoint_textbox.stored_string[] = "0"
            setpoint_textbox.displayed_string[] = "0"
        # elseif value != 1
        else
            multiplier = max(floor(Int, round(log(10, value - .1) - 1, digits = 8)), -1)

            setpoint_textbox.stored_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
            setpoint_textbox.displayed_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
        end
    else
        setpoint_textbox.stored_string[] = "5"
        setpoint_textbox.displayed_string[] = "5"
    end
end

on(setpoint_textbox.stored_string) do value
    setpoint.args[1][] = parse(Float64, value)
end

display(fig)      
