using GLMakie
using GLMakie.GLFW
using GLMakie: to_native

# dark_gray  = RGBf(0.5, 0.5, 0.5)
# mid_gray   = RGBf(0.7, 0.7, 0.7)
# light_gray = RGBf(0.85, 0.85, 0.85)
# off_white  = RGBf(0.9, 0.9, 0.9)

dark_gray  = RGBf(0.08, 0.08, 0.08)
mid_gray   = RGBf(0.16, 0.16, 0.16)
light_gray = RGBf(0.24, 0.24, 0.24)
off_white  = RGBf(0.5, 0.5, 0.5)

const table_textbox_attributes = Dict(
    :bordercolor               => light_gray,
    :bordercolor_focused       => light_gray,
    :bordercolor_hover         => light_gray,
    :borderwidth               => 1,
    :boxcolor                  => dark_gray,
    :boxcolor_focused          => mid_gray,
    :boxcolor_focused_invalid  => RGBf(1, 0.9, 0.9),
    :boxcolor_hover            => dark_gray,
    :cornerradius              => 0,
    :halign                    => :right,
    :height                    => 24,
    :placeholder               => " ",
    :textpadding               => (10, 8, 8, 4),
    :width                     => nothing,
)

const AXIS_ATTRS         = Dict(
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

const POP_UP_BUTTON_ATTRS       = Dict(
    :buttoncolor         => light_gray,
    :buttoncolor_active  => dark_gray,
    :buttoncolor_hover   => mid_gray,
    :cornerradius        => 4,
    :fontsize            => 16,
    :height              => 24,
    :labelcolor_hover    => :white,
    :width               => 96
)

# set_theme!(backgroundcolor = :gray95)

set_theme!(backgroundcolor = :gray12, textcolor = :gray80, fonts = (; 
    regular = "Times New Roman", 
    bold    = "Times New Roman Bold", 
    italic  = "Times New Roman Italic")
)

# set_theme!(backgroundcolor = :gray12, textcolor = :gray85, theme_latexfonts())

println(Makie.theme(:fonts))

# set_theme!(theme_dark())

fig_2 = Figure(size = (800, 400))

axis = Axis(fig_2[1, 1:5]; merge(AXIS_ATTRS, Dict(
    :title  => "Protocol",
    :xlabel => "Time [s]",
    :ylabel => "Setpoint"
))...)

# textbox = Textbox(fig_2[2, 1])
# add_button    = Button(fig_2[2, 2]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "ADD"))...)
# delete_button = Button(fig_2[2, 3]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "DELETE"))...)
# close_button  = Button(fig_2[2, 4]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "CLOSE"))...)

activate_button = Button(fig_2[2, 1]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Activate"))...)
import_button   = Button(fig_2[2, 2]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Import"))...)
export_button   = Button(fig_2[2, 3]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Export"))...)
clear_button    = Button(fig_2[2, 4]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Clear"))...)
close_button    = Button(fig_2[2, 5]; merge(POP_UP_BUTTON_ATTRS, Dict(:label => "Close"))...)


fig_2[3, 1:5] = table = GridLayout()

elements = []

for i in 1:10
    push!(elements, Textbox(table[1, i]; table_textbox_attributes...))
    push!(elements, Textbox(table[2, i]; table_textbox_attributes...))
end

colgap!(table, 0)
rowgap!(table, 0)

println(elements)

# elements[2].displayed_string[] = "test"


# Textbox(table[1, 1]; table_textbox_attributes...)
# Textbox(table[2, 1]; table_textbox_attributes...)
# rowgap!(table, 0)
# colgap!(table, 0)

# global event_number = 1
# global row_number = 1

# on(add.clicks) do _
#     global event_number = event_number + 1

#     Textbox(table[row_number, event_number]; table_textbox_attributes...)
#     Textbox(table[row_number + 1, event_number]; table_textbox_attributes...)

#     colgap!(table, 0)
#     rowgap!(table, 0)

#     for row in 1:((row_number-1) ÷ 2)
#         rowgap!(table, 2*row, 10)
#     end

#     if event_number == 10
#         # rowgap!(table, 2, 20)
#         global event_number = 0
#         global row_number = row_number + 2
#         # rowgap!(table, 0)
#         # rowgap!(table, row_number, 20)
#     end

#     println(fieldnames(typeof(fig_2)))
#     println(fig_2.content)
# end

# on(delete.clicks) do _
#     # pop!(fig_2.content)
#     # pop!(fig_2.content)
#     println(fig_2.content)
#     deleteat!(fig_2.content, 8)
#     println(fig_2.content)
# end

on(close_button.clicks) do _
    # if add.label[] == "Add"
    #     add.label[] = "Confirm"
    # else
    #     add.label[] = "Add"
    # end
    # fig_2.content = Any[]

    GLFW.SetWindowShouldClose(pop_up_window, true)
end

# global pop_up_screen = display(GLMakie.Screen(), fig_2)
# global pop_up_window = to_native(pop_up_screen)

fig_2
