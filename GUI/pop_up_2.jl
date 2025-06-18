using GLMakie
using GLMakie.GLFW
using GLMakie: to_native

include("functions.jl")

dark_gray  = RGBf(0.08, 0.08, 0.08)
mid_gray   = RGBf(0.16, 0.16, 0.16)
light_gray = RGBf(0.24, 0.24, 0.24)
off_white  = RGBf(0.5, 0.5, 0.5)

function int_range_validator(a::Integer, b::Integer)
    """
    Returns a validator function ensuring parsed integer s is in [a, b].
    """
    @assert a ≤ b

    return function(s::String)
        s_trim = strip(s)
        x = tryparse(Int, s_trim)
        return isempty(s_trim) || (x !== nothing && a ≤ x ≤ b)
    end
end

function range_validator(a::Integer, b::Integer)
    """
    Returns a validator function ensuring parsed integer s is in [a, b].
    """
    @assert a ≤ b

    return function(s::String)
        s_trim = strip(s)
        x = tryparse(Float64, s_trim)
        return isempty(s_trim) || (x !== nothing && a ≤ x ≤ b)
    end
end


TABLE_ATTRS = Dict(
    :bordercolor               => light_gray,
    :bordercolor_focused       => light_gray,
    :bordercolor_focused_invalid => light_gray,
    :bordercolor_hover         => light_gray,
    :borderwidth               => 1,
    :boxcolor                  => dark_gray,
    :boxcolor_focused          => mid_gray,
    :boxcolor_focused_invalid  => RGBf(0.22, 0.11, 0.11),
    :boxcolor_hover            => :gray12,
    :cornerradius              => 0,
    :halign                    => :right,
    :height                    => 24,
    :placeholder               => " ",
    :textcolor_placeholder     => :gray25,
    :textpadding               => (10, 8, 8, 4),
    :width                     => nothing,
)

AXIS_ATTRS = Dict(
    :alignmode        => Outside(10),
    :backgroundcolor  => dark_gray,
    :bottomspinecolor => light_gray,
    :halign           => :right,
    :height           => nothing,
    :leftspinecolor   => light_gray,
    :rightspinecolor  => light_gray,
    :spinewidth       => 0.5,
    :topspinecolor    => light_gray,
    :xgridcolor       => off_white,
    :xgridwidth       => 0.75,
    :xtickcolor       => off_white,
    :ygridcolor       => off_white,
    :ygridwidth       => 0.75,
    :ytickcolor       => off_white,
    :titlesize        => 16,
    :width            => nothing
)

TEXTBOX_ATTRS      = Dict(
    :bordercolor               => light_gray,
    :bordercolor_focused       => light_gray,
    :bordercolor_hover         => light_gray,
    :borderwidth               => 0.5,
    :boxcolor                  => dark_gray,
    :boxcolor_focused          => mid_gray,
    :boxcolor_focused_invalid  => RGBf(0.22, 0.11, 0.11),
    :boxcolor_hover            => light_gray,
    :cornerradius              => 2,
    :halign                    => :right,
    :height                    => 24,
    :placeholder               => " ",
    :textpadding               => (10, 8, 8, 4),
    :width                     => 48,
)

POP_UP_BUTTON_ATTRS = Dict(
    :buttoncolor         => light_gray,
    :buttoncolor_active  => dark_gray,
    :buttoncolor_hover   => mid_gray,
    :cornerradius        => 4,
    :fontsize            => 16,
    :height              => 24,
    :labelcolor_hover    => :white,
    :width               => 96
)

set_theme!(backgroundcolor = :gray12, textcolor = :gray80, fonts = (; 
    regular = "Times New Roman", 
    bold    = "Times New Roman Bold", 
    italic  = "Times New Roman Italic")
)

fig_3 = Figure(size = (300, 300))

Label(fig_3[1, 1], "Blue light & shutter controls", fontsize = 16)

label_grid = GridLayout(fig_3[2, 1])
buttons_grid = GridLayout(fig_3[3, 1])

Label(label_grid[1, 1], "Blue light duration [ms]", fontsize = 14, halign = :left)
duration = Textbox(label_grid[1, 2]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "100", :stored_string => "100"))...)
Label(label_grid[2, 1], "Delay [ms]", fontsize = 14, halign = :left)
dead_time = Textbox(label_grid[2, 2]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "50", :stored_string => "50"))...)
Label(label_grid[3, 1], "Frequency [Hz]", fontsize = 14, halign = :left)
freq = Textbox(label_grid[3, 2]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "1", :stored_string => "1"))...)

start_blue_button = Button(buttons_grid[1, 1]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Start"))...)
clear_blue_button = Button(buttons_grid[1, 2]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Clear"))...)

# global pop_up_screen = display(GLMakie.Screen(), fig_3)
# global pop_up_window = to_native(pop_up_screen)

on(start_blue_button.clicks) do _
    if connected[]
        command = "C 1\n"
        # command *= "S 1 S 1\n"
        # command *= "A 1 DP 1 $(parse(Float64, freq.stored_string[])) $(parse(Float64, duration.stored_string[])*parse(Float64, freq.stored_string[])/10)\n"
        # command *= "A 1 WT $(1000/parse(Float64, freq.stored_string[]) - parse(Float64, dead_time.stored_string[]))\n"
        # command *= "A 1 DP 3 $(parse(Float64, freq.stored_string[])) $((parse(Float64, duration.stored_string[]) + 2*parse(Float64, dead_time.stored_string[]))/parse(Float64, freq.stored_string[]))\n"
        # command *= "A 1 WT $(parse(Float64, duration.stored_string[]) + parse(Float64, dead_time.stored_string[]))\n"
        # command *= "A 1 DP 2 $(parse(Float64, freq.stored_string[])) $(100 - (parse(Float64, duration.stored_string[]) + 2*parse(Float64, dead_time.stored_string[]))*parse(Float64, freq.stored_string[])/10)\n"
        
        command *= "S 1 S 0\n"
        command *= "A 1 AO 1 0\n"
        command *= "A 1 WT $(parse(Float64, dead_time.stored_string[])*1000)\n"
        command *= "A 1 AO 2 5000\n"
        command *= "A 1 WT $(parse(Float64, duration.stored_string[])*1000)\n"
        command *= "A 1 AO 2 0\n"
        command *= "A 1 WT $(parse(Float64, dead_time.stored_string[])*1000)\n"
        command *= "A 1 AO 1 5000\n"
        command *= "A 1 WT $((1000/parse(Float64, freq.stored_string[]) - parse(Float64, duration.stored_string[]) - 2*parse(Float64, dead_time.stored_string[]))*1000)\n"

        
        command *= "A 0 TG 1\n"

        println(command)
        send_command(serial, command)
    else
        println("Not connected")
    end
end

on(clear_blue_button.clicks) do _
    if connected[]
        command = "C 1\nA 0 AS 1 0\nA 0 AS 2 0\n"

        println(command)
        send_command(serial, command)
    else
        println("Not connected")
    end
end

# fig_3
