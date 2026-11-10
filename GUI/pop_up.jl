using GLMakie
using GLMakie.GLFW
using GLMakie: to_native

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

function update_graph()
    empty!(axis)
    global protocol = []

    t_0 = parse(Int, delay.stored_string[])
    nb_repeats = parse(Int, repeats.stored_string[])

    cumtime = t_0

    vlines!(axis, 0, color = :transparent)
    vlines!(axis, cumtime, color=Makie.wong_colors()[6])

    # for _ in range(1, 2)
    for column in cell_values
        time, value = column

        if isnan(time)
            time = 0
        end

        lines!(axis, [cumtime, cumtime + time], [value, value], color=Makie.wong_colors()[1])
        if !isnan(value)
            vspan!(axis, cumtime, cumtime + time, color=Makie.wong_colors()[1], alpha=0.1)
        end
        cumtime += time
    end

    vlines!(axis, cumtime, color=Makie.wong_colors()[6])
    # end

    protocol = fill(NaN, round(Int, (cumtime - t_0)*nb_repeats + t_0+1))

    cumtime = t_0

    for _ in range(1, nb_repeats)
        for column in cell_values
            time, value = column

            if isnan(time)
                time = 0
            end

            protocol[round(Int, cumtime)+1:end] .= value
            cumtime += time
        end
    end

    # println(protocol)
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

TEXTBOX_ATTRS = Dict(
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

fig_2 = Figure(size = (600, 400))

fig_2[1, 1] = plotaxis = GridLayout()
fig_2[2, 1] = buttonaxis = GridLayout(tellwidth = false)
fig_2[3, 1] = table = GridLayout(tellwidth = false)

global axis = Axis(plotaxis[1, 1:3]; merge(AXIS_ATTRS, Dict(
    :alignmode => Inside(),
    :title     => "Protocol",
    :xlabel    => "Time [s]",
    :ylabel    => "Setpoint",
    :width     => nothing
))...)

# Box(buttonaxis[1:2, 1:3])

# activate_button = Button(fig_2[2, 1]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Activate"))...)
Label(buttonaxis[1, 1], "Repeats")
repeats = Textbox(buttonaxis[1, 2]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "1", :stored_string => "1", :validator => range_validator(0, 1000)))...)
Label(buttonaxis[2, 1], "Delay")
delay = Textbox(buttonaxis[2, 2]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "0", :stored_string => "0", :validator => range_validator(0, 3600)))...)

# Label(buttonaxis[1, 3], "Duty-cycle")
# duty_cycle = Textbox(buttonaxis[1, 4]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "0", :stored_string => "0", :validator => range_validator(0, 3600)))...)
# Label(buttonaxis[2, 3], "Dead-time")
# dead_time = Textbox(buttonaxis[2, 4]; merge(TEXTBOX_ATTRS, Dict(:displayed_string => "1", :stored_string => "1", :validator => range_validator(0, 1000)))...)
import_button   = Button(buttonaxis[1, 3]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Import"))...)
export_button   = Button(buttonaxis[2, 3]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Export"))...)
# clear_protocol_button    = Button(fig_2[2, 4]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Clear"))...)
clear_protocol_button    = Button(buttonaxis[1, 4]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Clear"))...)
active_button   = Button(buttonaxis[2, 4]; merge(POP_UP_BUTTON_ATTRS, Dict(
    :label             => "Activate",
    :buttoncolor       => mid_gray,
    :buttoncolor_hover => dark_gray
))...)

# fig_2[4, 1:3] = table = GridLayout()

global table_cell = []
global cell_values = []
global protocol = []

for i in 1:10
    push!(table_cell, Textbox(table[1, i]; merge(TABLE_ATTRS, Dict(:validator => int_range_validator(0, 3600)))...))
    push!(table_cell, Textbox(table[2, i]; merge(TABLE_ATTRS, Dict(:validator => range_validator(0, 10)))...))
    push!(cell_values, [NaN, NaN])
end

colgap!(table, 0)
rowgap!(table, 0)

# println(cell_values)

on(repeats.stored_string) do _
    update_graph()
end

on(delay.stored_string) do _
    update_graph()
end

for (i, textbox) in enumerate(table_cell)
    on(textbox.stored_string) do val
        if val == " "
            cell_values[cld(i, 2)][mod1(i, 2)] = NaN
        else
            cell_values[cld(i, 2)][mod1(i, 2)] = parse(Float64, val)
        end

        update_graph()
    end
end

# on(clear_protocol_button.clicks) do _
#     for textbox in table_cell
#         textbox.stored_string[] = " "
#     end
# end

on(clear_protocol_button.clicks) do _
    # GLFW.SetWindowShouldClose(pop_up_window, true)
    
    while !isempty(table_cell)
        delete!(table_cell[end])
        pop!(table_cell)
    end

    while !isempty(cell_values)
        pop!(cell_values)
    end

    global cell_values = []

    for i in 1:10
        push!(table_cell, Textbox(table[1, i]; merge(TABLE_ATTRS, Dict(:validator => int_range_validator(0, 3600)))...))
        push!(table_cell, Textbox(table[2, i]; merge(TABLE_ATTRS, Dict(:validator => range_validator(0, 10)))...))
        push!(cell_values, [NaN, NaN])
    end

    for (i, textbox) in enumerate(table_cell)
        on(textbox.stored_string) do val
            if val == " "
                cell_values[cld(i, 2)][mod1(i, 2)] = NaN
            else
                cell_values[cld(i, 2)][mod1(i, 2)] = parse(Float64, val)
            end

            update_graph()
        end
    end

    update_graph()
end

on(active_button.clicks) do _
    if protocol_active[]
        active_button.buttoncolor[] = mid_gray
        active_button.label[] = "Activate"
        protocol_active[] = false
    else
        active_button.buttoncolor[] = light_gray
        active_button.label[] = "Deactivate"
        protocol_active[] = true
    end
end

# global pop_up_screen = display(GLMakie.Screen(), fig_2)
# global pop_up_window = to_native(pop_up_screen)

# fig_2
