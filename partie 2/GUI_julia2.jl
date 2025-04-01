using LibSerialPort
using GLMakie
using Observables

try
    LibSerialPort.close(serial)
catch
    print()
end

function pwm_vector(periods, resoltion, duty_cycle, delay, intensity)
    vector = zeros(resoltion)
    num_ones = round(Int, resoltion * duty_cycle / 100/periods)
    start_index = round(Int, resoltion * delay / 100/periods)

    for period in 1:periods
        for i in 0:num_ones-1
            vector[round(Int, (start_index + i + (period + 1)/periods*resolution)) % resoltion + 1] = intensity
        end
    end

    return vector
end

function camera_vector(periods, resolution, duty_cycle, delay, images, duration)
    vector = zeros(resolution)

    num_ones = round(Int, duration*resolution * duty_cycle / 5000/periods)
    num_zeros = round(Int, (100 - duration)*resolution * duty_cycle / 5000/periods)
    start_index = round(Int, resolution * delay / 100/periods)
    pulse_duration = num_ones/images
    pulse_interval = num_zeros/images

    for period in 1:periods
        for i in 0:images-1
            for j in 0:div(pulse_duration, 2) - 1
                index = round(Int, (start_index + pulse_interval/4 + i*(pulse_duration+pulse_interval)/2 + j + (period + 1)/periods*resolution)) % resolution + 1

                if vector[index] == 100
                    vector[index + 1] = 100
                else
                    vector[index] = 100
                end
            end
        end
    end

    return vector
end

function add_to_monitor(string)
    initial_string = monitor.text[]
    lines = split(initial_string, "\n")
    push!(lines, split(string, "\r\n")[1])
    lines = filter(!isempty, lines)

    nb_lines = 17

    if length(lines) > nb_lines
        monitor.text[] = join(lines[end - (nb_lines-1):end], "\n")
    else
        monitor.text[] = join(lines, "\n")
    end
end

function connect_to_port(portname)
    try
        # if isopen("/dev/tty."*portname)
        #     close("/dev/tty."*portname)
        # end
        serial = LibSerialPort.open("/dev/tty."*portname, 115200)
        LibSerialPort.set_read_timeout(serial, 3)
        sleep(0.1)
        add_to_monitor("Connection successful")
        println("Connection successful")
        connected[] = true
        connect.label[] = "Disconnect"

        return serial
    catch e
        # throw(e)
        add_to_monitor("Connection error")
        println("Connection error")
    end
end

function send_command(serial, command)
    println("sending signal")
    write(serial, command)
    sleep(0.1)
    response = strip(read(serial, String))
    println(string(response))

    add_to_monitor("$response")

    return response
end

function create_LED_command(frequency, channel_1, duty_cycle_1, delay_1, intensity_1, channel_2, duty_cycle_2, delay_2, intensity_2)
    command = ""

    if channel_1 || channel_2
        channel_1_command = "A 1 AP 1 $frequency $duty_cycle_1 0 $(intensity_1*100)\n"
        channel_2_command = "A 1 AP 2 $frequency $duty_cycle_2 0 $(intensity_2*100)\n"
        # channel_1_command = "A 1 AO 1 $(intensity_1*100)\nA 1 DP 1 $frequency $duty_cycle_1\n"
        # channel_2_command = "A 1 AO 2 $(intensity_2*100)\nA 1 DP 2 $frequency $duty_cycle_2\n"
    
        if delay_1 < delay_2
            command *= "A 1 WT $(round(Int, delay_1/frequency*10000))\n"

            if channel_1
                command *= channel_1_command
            end

            command *= "A 1 WT $(round(Int, (delay_2 - delay_1)/frequency*10000))\n"
            
            if channel_2
                command *= channel_2_command
            end

            command *= "A 1 WT $(round(Int, (100 - delay_2)/frequency*10000))\n"
        else
            command *= "A 1 WT $(round(Int, delay_1/frequency*10000))\n"
            
            if channel_2
                command *= channel_2_command
            end

            command *= "A 1 WT $(round(Int, (delay_1 - delay_2)/frequency*10000))\n"

            if channel_1
                command *= channel_1_command
            end

            command *= "A 1 WT $(round(Int, (100 - delay_1)/frequency*10000))\n"
        end
    else
        println("No command available")
    end

    return command
end

function create_camera_command(frequency, channel_1, duty_cycle_1, delay_1, images_1, duration_1, channel_2, duty_cycle_2, delay_2, images_2, duration_2)
    command = ""

    if channel_1 || channel_2
        channel_1_command = ""
        channel_2_command = ""

        for _ in 1:images_1
            channel_1_command *= "A 1 WT $(round(Int, (100 - duration_1)*duty_cycle_1/images_1/frequency*50))\n"
            channel_1_command *= "A 1 DP 1 $frequency $(duration_1*duty_cycle_1/images_1/100)\n"
            channel_1_command *= "A 1 WT $(round(Int, (100 + duration_1)*duty_cycle_1/images_1/frequency*50))\n"
            # channel_1_command *= "A 1 WT $(round(Int, (duty_cycle_1/images_1/4)/frequency*10000))\n"
            # channel_1_command *= "A 1 DP 1 $frequency $(duty_cycle_1/images_1/2)\n"
            # channel_1_command *= "A 1 WT $(round(Int, (3*duty_cycle_1/images_1/4)/frequency*10000))\n"
        end

        for _ in 1:images_2
            channel_2_command *= "A 1 WT $(round(Int, (100 - duration_2)*duty_cycle_2/images_2/frequency*50))\n"
            channel_2_command *= "A 1 DP 2 $frequency $(duration_2*duty_cycle_2/images_2/100)\n"
            channel_2_command *= "A 1 WT $(round(Int, (100 + duration_2)*duty_cycle_2/images_2/frequency*50))\n"
            # channel_2_command *= "A 1 WT $(round(Int, (duty_cycle_2/images_2/4)/frequency*10000))\n"
            # channel_2_command *= "A 1 DP 2 $frequency $(duty_cycle_2/images_2/2)\n"
            # channel_2_command *= "A 1 WT $(round(Int, (3*duty_cycle_2/images_2/4)/frequency*10000))\n"
        end
    
        if delay_1 < delay_2
            command *= "A 1 WT $(round(Int, delay_1/frequency*10000))\n"

            if channel_1
                command *= channel_1_command
            end

            command *= "A 1 WT $(round(Int, (delay_2 - delay_1 - duty_cycle_1)/frequency*10000))\n"
            
            if channel_2
                command *= channel_2_command
            end

            command *= "A 1 WT $(round(Int, (100 - delay_2 - duty_cycle_2)/frequency*10000))\n"
        else
            command *= "A 1 WT $(round(Int, delay_1/frequency*10000))\n"
            
            if channel_2
                command *= channel_2_command
            end

            command *= "A 1 WT $(round(Int, (delay_1 - delay_2 - duty_cycle_2)/frequency*10000))\n"

            if channel_1
                command *= channel_1_command
            end

            command *= "A 1 WT $(round(Int, (100 - delay_1 - duty_cycle_1)/frequency*10000))\n"
        end
    else
        println("No command available")
    end

    return command
end

function create_command(frequency, channel_1, duty_cycle_1, delay_1, intensity_1, images_1, duration_1, channel_2, duty_cycle_2, delay_2, intensity_2, images_2, duration_2)
    # LED_command = create_LED_command(frequency, channel_1, duty_cycle_1, delay_1, intensity_1, channel_2, duty_cycle_2, delay_2, intensity_2)
    LED_command = ""
    camera_command = create_camera_command(frequency, channel_1, duty_cycle_1, delay_1, images_1, duration_1, channel_2, duty_cycle_2, delay_2, images_2, duration_2)
    # camera_command = ""

    command = "C Z\n"
    command *= "C R\n"
    command *= "S 1 S 1\n"
    command *= LED_command
    command *= camera_command

    command *= "A 0 TG 1\n"

    return command
end

set_theme!(backgroundcolor = :gray95)

fig = Figure(size = (700, 800))

fig[1, 1] = plotgrid = GridLayout(tellwidth=true)
fig[2, 1] = blockgrid = GridLayout(tellwidth=false)
blockgrid[1, 1] = buttongrid = GridLayout(tellwidth=false)
blockgrid[1, 2] = slidergrid = GridLayout(tellwidth=false)
buttongrid[2, 1] = frequencygrid = GridLayout(tellwidth=false)
buttongrid[2, 2] = fpsgrid = GridLayout(tellwidth=false)

colsize!(blockgrid, 1, 200)
colgap!(blockgrid, 1, 64)

frequency_label = Label(buttongrid[1, 1], text="Frequency [Hz] :", fontsize = 16, justification=:right)

frequency_textbox = Textbox(frequencygrid[1:2, 1], placeholder = " ", validator = Float64, width = 96, height = 32, cornerradius = 0, borderwidth=0,
    boxcolor=:white, boxcolor_focused=:white, boxcolor_focused_invalid=RGBf(1, 0.9, 0.9), boxcolor_hover=RGBf(0.975, 0.975, 0.975))

frequency_up = Button(frequencygrid[1, 2], label="▲", font=:bold, fontsize=8, width = 32, height = 17, cornerradius = 0, valign=:bottom, strokewidth=1,
    buttoncolor=:gray95, buttoncolor_hover=RGBf(0.8, 0.8, 0.8), buttoncolor_active=RGBf(0.4, 0.4, 0.4), labelcolor_active=:black, strokecolor=RGBf(0.4, 0.4, 0.4))

frequency_down = Button(frequencygrid[2, 2], label="▼", font=:bold, fontsize=8, width = 32, height = 17, cornerradius = 0, valign=:top, strokewidth=1,
    buttoncolor=:gray95, buttoncolor_hover=RGBf(0.8, 0.8, 0.8), buttoncolor_active=RGBf(0.4, 0.4, 0.4), labelcolor_active=:black, strokecolor=RGBf(0.4, 0.4, 0.4))

Box(frequencygrid[1:2, 1:2], color=:transparent, width = 128, height = 34, cornerradius = 0)

frequency_textbox.displayed_string[] = "10.0"
frequency_textbox.stored_string[] = "10.0"

fps_label = Label(buttongrid[1, 2], text="Max fps:", fontsize = 16, justification=:right)

fps_textbox = Textbox(fpsgrid[1:2, 1], placeholder = " ", validator = UInt16, width = 96, height = 32, cornerradius = 0, borderwidth=0,
    boxcolor=:white, boxcolor_focused=:white, boxcolor_focused_invalid=RGBf(1, 0.9, 0.9), boxcolor_hover=RGBf(0.975, 0.975, 0.975))

fps_up = Button(fpsgrid[1, 2], label="▲", font=:bold, fontsize=8, width = 32, height = 17, cornerradius = 0, valign=:bottom, strokewidth=1,
    buttoncolor=:gray95, buttoncolor_hover=RGBf(0.8, 0.8, 0.8), buttoncolor_active=RGBf(0.4, 0.4, 0.4), labelcolor_active=:black, strokecolor=RGBf(0.4, 0.4, 0.4))

fps_down = Button(fpsgrid[2, 2], label="▼", font=:bold, fontsize=8, width = 32, height = 17, cornerradius = 0, valign=:top, strokewidth=1,
    buttoncolor=:gray95, buttoncolor_hover=RGBf(0.8, 0.8, 0.8), buttoncolor_active=RGBf(0.4, 0.4, 0.4), labelcolor_active=:black, strokecolor=RGBf(0.4, 0.4, 0.4))

Box(fpsgrid[1:2, 1:2], color=:transparent, width = 128, height = 34, cornerradius = 0)

fps_textbox.displayed_string[] = "1000"
fps_textbox.stored_string[] = "1000"

colgap!(frequencygrid, 1, Fixed(-1))
rowgap!(frequencygrid, 1, Fixed(0))
colgap!(fpsgrid, 1, Fixed(-1))
rowgap!(fpsgrid, 1, Fixed(0))
colgap!(buttongrid, 1, Fixed(18))

start = Button(buttongrid[3, 1], label="START", font=:bold, fontsize=16, width = 128, height = 32, cornerradius = 0, valign=:bottom, strokewidth=2,
    buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), strokecolor = :black, labelcolor_active=:black)
    
stop = Button(buttongrid[3, 2], label="STOP", font=:bold, fontsize=16, width = 128, height = 32, cornerradius = 0, valign=:bottom, strokewidth=2,
    buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), strokecolor = :black, labelcolor_active=:black)

portname_textbox = Textbox(buttongrid[4, 1], placeholder = "Port name...", width = 128, height = 32, cornerradius = 0, borderwidth=1,
    boxcolor=:white, boxcolor_focused=:white, boxcolor_focused_invalid=RGBf(1, 0.9, 0.9), boxcolor_hover=RGBf(0.975, 0.975, 0.975), bordercolor = :black, bordercolor_focused = :black, bordercolor_hover = :black)

portname_textbox.displayed_string[] = "usbmodem102"
portname_textbox.stored_string[] = "usbmodem102"

connect = Button(buttongrid[4, 2], label="Connect", font=:bold, fontsize=16, width = 128, height = 32, cornerradius = 0, valign=:bottom, strokewidth=2,
    buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), strokecolor = :black, labelcolor_active=:black)

monitor_box = Box(buttongrid[5, 1:2], color=:white, width = 274, height = 296,cornerradius = 0)

global monitor = Label(buttongrid[5, 1:2], "\n", width = 274, height = 120, justification = :left, valign = :bottom, padding = (8, 8, 8, 8), tellheight = false, tellwidth = false)


ch1_label = Label(slidergrid[1, 1][1, 1], text="Channel 1", font=:bold, fontsize = 18)

ch1_check = Toggle(slidergrid[1, 1][1, 2], active = true, rimfraction=0.24, 
    buttoncolor=RGBf(0.4, 0.4, 0.4), framecolor_active=RGBf(0.6, 0.6, 0.6), framecolor_inactive=RGBf(0.8, 0.8, 0.8))

ch1_menu = Menu(slidergrid[1, 1][1, 3], options = ["Square", "Tri", "Exponential"], width = 96, height = 34, default = "Square", selection_cell_color_inactive=:white)
Box(slidergrid[1, 1][1, 3], color=:transparent, width = 96, height = 34, cornerradius = 0)

ch1_sliders = SliderGrid(
    slidergrid[1, 1][2, :],
    (label = "Duty cycle", range = 0:1:100, format = "{:d} %", startvalue = 40, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Delay", range = 0:1:100, format = "{:d} %", startvalue = 5, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Intensity", range = 0:1:100, format = "{:d} %", startvalue = 100, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Images", range = 0:1:40, format = " {:d} (10000 fps)", startvalue = 5, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Pulse\nduration", range = 0:1:100, format = "{:d} % (999 ms)", startvalue = 50, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    width = 350)
ch1_sliders.layout.content[12].content.text = "5 (125 fps)"
ch1_sliders.layout.content[15].content.text = "50 % (4.0 ms)"

ch2_label = Label(slidergrid[3, 1][1, 1], text="Channel 2", font=:bold, fontsize = 18)

ch2_check = Toggle(slidergrid[3, 1][1, 2], active = true, rimfraction=0.24,
    buttoncolor=RGBf(0.4, 0.4, 0.4), framecolor_active=RGBf(0.6, 0.6, 0.6), framecolor_inactive=RGBf(0.8, 0.8, 0.8))

ch2_sliders = SliderGrid(
    slidergrid[3, 1][2, :],
    (label = "Duty cycle", range = 0:1:100, format = "{:d} %", startvalue = 40, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Delay", range = 0:1:100, format = "{:d} %", startvalue = 55, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Intensity", range = 0:1:100, format = "{:d} %", startvalue = 100, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Images", range = 0:1:40, format = " {:d} (10000 fps)", startvalue = 2, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    (label = "Pulse\nduration", range = 0:1:100, format = "{:d} % (999 ms)", startvalue = 50, linewidth = 14,
        color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false),
    width = 350)
ch2_sliders.layout.content[12].content.text = "2 (50 fps)"
ch2_sliders.layout.content[15].content.text = "50 % (10.0 ms)"

rowsize!(slidergrid, 2, Fixed(20))

periods = 2
resolution = 20000
t = range(0, periods, resolution)
low = zeros(resolution)
LED_1 = Observable(zeros(resolution))
LED_1[] = pwm_vector(periods, resolution, 40, 5, 100)
LED_2 = Observable(zeros(resolution))
LED_2[] = pwm_vector(periods, resolution, 40, 55, 100)
camera_1 = Observable(zeros(resolution))
camera_1[] = camera_vector(periods, resolution, 40, 5, 5, 50)
camera_2 = Observable(zeros(resolution))
camera_2[] = camera_vector(periods, resolution, 40, 55, 2 ,50)
error = Observable(zeros(resolution))
error[] = camera_1[].*camera_2[]./100
global serial = nothing
global connected = Observable(false)

# global axis = plotgrid[1, 1] = Axis(fig, title="Period", xlabel = "Time", ylabel = "Intensity [%]", titlesize=16, aspect = 4, limits = (-0.1, periods + 0.1, -10, 110), 
#     xticks = 0:periods, xtickformat = "{:d} T")
global axis = plotgrid[1, 1] = Axis(fig, title="Signal", ylabel = "Intensity [%]", titlesize=16, aspect = 4, limits = (-0.1, periods + 0.1, -10, 110), 
    xticks = 0:periods, xtickformat = "{:d} T")

band!(axis, t, LED_1, low, color=Makie.wong_colors()[1], alpha = 0.25)
band!(axis, t, LED_2, low, color=Makie.wong_colors()[2], alpha = 0.25)
camera_1_band = band!(axis, t, camera_1, low, color=:black, alpha = 0.5)
camera_2_band = band!(axis, t, camera_2, low, color=:black, alpha = 0.5)
lines!(axis, t, LED_1, color=Makie.wong_colors()[1])
lines!(axis, t, LED_2, color=Makie.wong_colors()[2])
error_band = band!(axis, t, error, low, color=:red2, alpha = 0.5)
lines!(axis, t, low, color=:black)

on(frequency_textbox.stored_string) do string
    if parse(Float64, string) < 0.04
        frequency_textbox.displayed_string[] = "0.04"
        frequency_textbox.stored_string[] = "0.04"
    elseif parse(Float64, string) > 1000000
        frequency_textbox.displayed_string[] = "1000000.0"
        frequency_textbox.stored_string[] = "1000000.0"
    end

    max_images_1 = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, string)*ch1_sliders.sliders[1].value[]/100)
    max_images_2 = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, string)*ch2_sliders.sliders[1].value[]/100)
    ch1_sliders.sliders[4].range = 0:1:max_images_1
    ch2_sliders.sliders[4].range = 0:1:max_images_2
end

on(frequency_up.clicks) do _
    original_frequency = parse(Float64, frequency_textbox.stored_string[])

    # if log10(original_frequency) % 1 == 0.1
    #     frequency = ceil(Int, original_frequency/(10^floor(Int, log10(original_frequency) - 1)) + 1e-9)*10^floor(Int, log10(original_frequency) - 1)
    # else
    #     frequency = ceil(Int, original_frequency/(10^floor(Int, log10(original_frequency))) + 1e-9)*10^floor(Int, log10(original_frequency))
    # end

    frequency = ceil(Int, original_frequency/(10^floor(Int, log10(original_frequency))) + 1e-9)*10^floor(Int, log10(original_frequency))

    if frequency < 0.04
        frequency_textbox.displayed_string[] = "0.04"
        frequency_textbox.stored_string[] = "0.04"
    elseif frequency > 1000000
        frequency_textbox.displayed_string[] = "1000000.0"
        frequency_textbox.stored_string[] = "1000000.0"
    else
        frequency_textbox.displayed_string[] = string(frequency)*".0"
        frequency_textbox.stored_string[] = string(frequency)*".0"
    end
end

on(frequency_down.clicks) do _
    # frequency = floor(Int, parse(Float64, frequency_textbox.displayed_string[]) - 1e-9)
    original_frequency = parse(Float64, frequency_textbox.stored_string[])

    if log10(original_frequency) % 1 == 0
        frequency = floor(Int, original_frequency/(10^floor(Int, log10(original_frequency) - 1)) - 1e-9)*10^floor(Int, log10(original_frequency) - 1)
    else
        frequency = floor(Int, original_frequency/(10^floor(Int, log10(original_frequency))) - 1e-9)*10^floor(Int, log10(original_frequency))
    end

    if frequency < 0.04
        frequency_textbox.displayed_string[] = "0.04"
        frequency_textbox.stored_string[] = "0.04"
    elseif frequency > 1000000
        frequency_textbox.displayed_string[] = "1000000.0"
        frequency_textbox.stored_string[] = "1000000.0"
    else
        frequency_textbox.displayed_string[] = string(frequency)*".0"
        frequency_textbox.stored_string[] = string(frequency)*".0"
    end
end

on(fps_textbox.stored_string) do string
    max_images_1 = floor(Int, parse(Int, string)/parse(Float64, frequency_textbox.stored_string[])*ch1_sliders.sliders[1].value[]/100)
    max_images_2 = floor(Int, parse(Int, string)/parse(Float64, frequency_textbox.stored_string[])*ch2_sliders.sliders[1].value[]/100)
    ch1_sliders.sliders[4].range = 0:1:max_images_1
    ch2_sliders.sliders[4].range = 0:1:max_images_2
end

on(fps_up.clicks) do _
    original_fps = parse(Float64, fps_textbox.stored_string[])

    fps = ceil(Int, original_fps/(10^floor(Int, log10(original_fps))) + 1e-9)*10^floor(Int, log10(original_fps))

    fps_textbox.displayed_string[] = string(fps)
    fps_textbox.stored_string[] = string(fps)
end

on(fps_down.clicks) do _
    original_fps = parse(Float64, fps_textbox.stored_string[])

    if log10(original_fps) % 1 == 0
        fps = floor(Int, original_fps/(10^floor(Int, log10(original_fps) - 1)) - 1e-9)*10^floor(Int, log10(original_fps) - 1)
    else
        fps = floor(Int, original_fps/(10^floor(Int, log10(original_fps))) - 1e-9)*10^floor(Int, log10(original_fps))
    end

    fps_textbox.displayed_string[] = string(fps)
    fps_textbox.stored_string[] = string(fps)
end

on(start.clicks) do _
    command = create_command(round(parse(Float64, frequency_textbox.displayed_string[]), digits = 2), 
            ch1_check.active[], ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[], 
            ch2_check.active[], ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
    println(command)

    if connected[]
        if error[] == zeros(resolution)
            command = create_command(round(parse(Float64, frequency_textbox.displayed_string[]), digits = 2), 
            ch1_check.active[], ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[], 
            ch2_check.active[], ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
            # println(command)

            send_command(serial, command)
        else
            add_to_monitor("Camera superposition error")
            println("Camera superposition error")
        end
    else
        add_to_monitor("Not connected")
        println("Not connected")
    end
end

on(stop.clicks) do _
    if connected[]
        send_command(serial, "C Z\n")
    else
        add_to_monitor("Not connected")
        println("Not connected")
    end
end

on(connect.clicks) do _
    if !connected[]
        global serial = connect_to_port(portname_textbox.stored_string[])
    else
        LibSerialPort.close(serial)
        add_to_monitor("Disconnection successful")
        println("Disconnection successful")
        connect.label[] = "Connect"
        connected[] = false
    end
end

on(ch1_sliders.sliders[1].value) do value
    if ch1_check.active[]
        LED_1[] = pwm_vector(periods, resolution, value, ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[])

        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/periods
            camera_1[] = pwm_vector(periods, resolution, value, ch1_sliders.sliders[2].value[], 100)
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(periods, resolution, value, ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[])
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    end

    max_images = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, frequency_textbox.stored_string[])*value/100)
    ch1_sliders.sliders[4].range = 0:1:max_images

    fps = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch1_sliders.sliders[2].value) do value
    if ch1_check.active[]
        LED_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], value, ch1_sliders.sliders[3].value[])

        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/periods
            camera_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], value, 100)
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(periods, resolution, ch1_sliders.sliders[1].value[], value, ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[])
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    end
end

on(ch1_sliders.sliders[3].value) do value
    if ch1_check.active[]
        LED_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], value)
        error[] = camera_1[].*camera_2[]./100
    end
end

on(ch1_sliders.sliders[4].value) do value
    if ch1_check.active[]
        if value/ch1_sliders.sliders[1].value[] >= 2/periods
            camera_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], 100)
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], value, ch1_sliders.sliders[5].value[])
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        
        error[] = camera_1[].*camera_2[]./100
    end

    ch1_sliders.layout.content[12].content.text = "$value ($(round(Int, 100*value*parse(Float64, frequency_textbox.stored_string[])/ch1_sliders.sliders[1].value[])) fps)"

    fps = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch1_sliders.sliders[5].value) do value
    if ch1_check.active[]
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/periods
            camera_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], 100)
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[], value)
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        
        error[] = camera_1[].*camera_2[]./100
    end

    fps = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch2_sliders.sliders[1].value) do value
    if ch2_check.active[]
        LED_2[] = pwm_vector(periods, resolution, value, ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[])

        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/periods
            camera_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 100)
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(periods, resolution, value, ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
            camera_2_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    end

    max_images = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, frequency_textbox.stored_string[])*value/100)
    ch2_sliders.sliders[4].range = 0:1:max_images

    fps = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch2_sliders.sliders[2].value) do value
    if ch2_check.active[]
        LED_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], value, ch2_sliders.sliders[3].value[])

        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/periods
            camera_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 100)
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(periods, resolution, ch2_sliders.sliders[1].value[], value, ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
            camera_2_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    end
end

on(ch2_sliders.sliders[3].value) do value
    if ch2_check.active[]
        LED_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], value)
        error[] = camera_1[].*camera_2[]./100
    end
end

on(ch2_sliders.sliders[4].value) do value
    if ch2_check.active[]
        if value/ch2_sliders.sliders[1].value[] >= 2/periods
            camera_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 100)
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], value, ch2_sliders.sliders[5].value[])
            camera_2_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    end

    ch2_sliders.layout.content[12].content.text = "$value ($(round(Int, 100*value*parse(Float64, frequency_textbox.stored_string[])/ch2_sliders.sliders[1].value[])) fps)"

    fps = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch2_sliders.sliders[5].value) do value
    if ch2_check.active[]
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/periods
            camera_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 100)
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[], value)
            camera_2_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        
        error[] = camera_1[].*camera_2[]./100
    end

    fps = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)

    if fps >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps)) μs)"
    end
end

on(ch1_check.active) do active
    if active
        LED_1[] = pwm_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[])

        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/periods
            camera_1[] = LED_1[]
            camera_1_band.alpha = 0.25
            error_band.alpha = 0.25
        else
            camera_1[] = camera_vector(periods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[])
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    else
        LED_1[] = zeros(resolution)
        camera_1[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100
    end
end

on(ch2_check.active) do active
    if active
        LED_2[] = pwm_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[])

        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/periods
            camera_2[] = LED_2[]
            camera_2_band.alpha = 0.25
        else
            camera_2[] = camera_vector(periods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
            camera_2_band.alpha = 0.5
        end

        error[] = camera_1[].*camera_2[]./100
    else
        LED_2[] = zeros(resolution)
        camera_2[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100
    end
end

fig
