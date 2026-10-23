########################################################################
# LIBRARY IMPORTS & INITIAL SETUP
########################################################################
using LibSerialPort
using GLMakie
using Observables

# Close any previously opened serial connection (if any)
try
    LibSerialPort.close(serial)
catch
    # Nothing to close
    println("No serial connection to close.")
end

########################################################################
# SIGNAL GENERATION FUNCTIONS
########################################################################

"""
    ch_vector_gen(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates a square signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function square_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity
        end
    end
    return signal
end

"""
    saw_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates a sawtooth signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function saw_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity*i/numOn
        end
    end
    return signal
end

"""
    rsaw_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates a reverse sawtooth signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function rsaw_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity*(numOn-i)/numOn
        end
    end
    return signal
end

"""
    tri_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates a triangle signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function tri_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity*(numOn-abs(numOn-2*i))/numOn
        end
    end
    return signal
end

"""
    rtri_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates a reverse triangle signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function rtri_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity*(abs(numOn-2*i))/numOn
        end
    end
    return signal
end

"""
    exp_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)

Generates an exponential signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Number of samples in one period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay (percent of period) before the pulse.
  - intensity: Intensity value for the pulses.

Returns a vector of length `resolution` with pulses of given intensity.
"""
function exp_vector(numPeriods, resolution, dutyPercent, delayPercent, intensity)
    signal = zeros(resolution)
    numOn = round(Int, resolution * dutyPercent / 100 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    
    for period in 1:numPeriods
        for i in 0:numOn-1
            # Compute the index and use modulo for wrapping
            idx = round(Int, (startIdx + i + (period + 1) / numPeriods * resolution)) % resolution + 1
            signal[idx] = intensity*exp(-3*i/numOn)
        end
    end
    return signal
end

"""
    camera_vector(numPeriods, resolution, dutyPercent, delayPercent, numImages, durationPercent)

Generates a camera trigger signal vector.

Arguments:
  - numPeriods: Number of periods.
  - resolution: Samples per period.
  - dutyPercent: Duty cycle in percent.
  - delayPercent: Delay in percent.
  - numImages: Number of trigger pulses per period.
  - durationPercent: Parameter affecting pulse duration.

Returns a vector of length `resolution` with camera trigger pulses.
"""
function camera_vector(numPeriods, resolution, dutyPercent, delayPercent, numImages, durationPercent)
    signal = zeros(resolution)
    
    numOn = round(Int, durationPercent * resolution * dutyPercent / 5000 / numPeriods)
    numOff = round(Int, (100 - durationPercent) * resolution * dutyPercent / 5000 / numPeriods)
    startIdx = round(Int, resolution * delayPercent / 100 / numPeriods)
    pulseDur = numOn / numImages
    pulseIntvl = numOff / numImages

    for period in 1:numPeriods
        for i in 0:numImages-1
            for j in 0:div(pulseDur, 2) - 1
                idx = round(Int, (startIdx + pulseIntvl/4 + i * (pulseDur + pulseIntvl)/2 + j + (period + 1)/numPeriods * resolution)) % resolution + 1
                # If already set to 100 then set the next index also to 100; otherwise, set current to 100.
                if signal[idx] == 120
                    signal[idx + 1] = 120
                else
                    signal[idx] = 120
                end
            end
        end
    end
    return signal
end

########################################################################
# COMMUNICATION & MONITOR FUNCTIONS
########################################################################

"""
    add_to_monitor(msg::String)

Appends a message to the GUI monitor (keeps only the last 17 lines).
"""
function add_to_monitor(msg::String)
    currentText = monitor.text[]
    lines = split(currentText, "\n")
    push!(lines, split(msg, "\r\n")[1])
    lines = filter(!isempty, lines)
    MAX_LINES = 17
    if length(lines) > MAX_LINES
        monitor.text[] = join(lines[end - (MAX_LINES - 1):end], "\n")
    else
        monitor.text[] = join(lines, "\n")
    end
end

"""
    connect_to_port(portName)

Attempts to open a serial connection on "/dev/tty.<portName>".
If successful, updates the monitor and GUI state.

Returns the serial object (or nothing on failure).
"""
function connect_to_port(portName)
    try
        ser = LibSerialPort.open("/dev/tty." * portName, 115200)
        LibSerialPort.set_read_timeout(ser, 3)
        sleep(0.1)
        add_to_monitor("Connection successful")
        println("Connection successful")
        connected[] = true
        connect.label[] = "Disconnect"
        return ser
    catch e
        add_to_monitor("Connection error")
        println("Connection error")
    end
end

"""
    send_command(serialConn, commandStr)

Sends the command string over the serial connection.
Waits briefly, reads the response, prints and appends it to the monitor.
"""
function send_command(serialConn, commandStr)
    println("sending signal")
    println(commandStr)
    write(serialConn, commandStr)
    sleep(0.1)
    response = strip(read(serialConn, String))
    println(response)
    add_to_monitor("$response")
    return response
end

########################################################################
# COMMAND CREATION FUNCTIONS
########################################################################

"""
    create_LED_command(freq, ch1Enabled, duty1, delay1, intensity1,
                       ch2Enabled, duty2, delay2, intensity2)

Builds the LED command string (if any channel is enabled).

Note: Some alternate command strings are commented.
"""
function create_LED_command(freq, ch1Enabled, duty1, delay1, intensity1,
                            ch2Enabled, duty2, delay2, intensity2)
    cmd = ""
    if ch1Enabled || ch2Enabled
        cmd_ch1 = "A 1 AP 1 $freq $duty1 0 $(intensity1 * 100)\n"
        cmd_ch2 = "A 1 AP 2 $freq $duty2 0 $(intensity2 * 100)\n"
        if delay1 < delay2
            cmd *= "A 1 WT $(round(Int, delay1/freq * 10000))\n"
            if ch1Enabled
                cmd *= cmd_ch1
            end
            cmd *= "A 1 WT $(round(Int, (delay2 - delay1)/freq * 10000))\n"
            if ch2Enabled
                cmd *= cmd_ch2
            end
            cmd *= "A 1 WT $(round(Int, (100 - delay2)/freq * 10000))\n"
        else
            cmd *= "A 1 WT $(round(Int, delay2/freq * 10000))\n"
            if ch2Enabled
                cmd *= cmd_ch2
            end
            cmd *= "A 1 WT $(round(Int, (delay1 - delay2)/freq * 10000))\n"
            if ch1Enabled
                cmd *= cmd_ch1
            end
            cmd *= "A 1 WT $(round(Int, (100 - delay1)/freq * 10000))\n"
        end
    else
        println("No command available")
    end
    return cmd
end

"""
    create_camera_command(freq, ch1Enabled, duty1, delay1, images1, duration1,
                          ch2Enabled, duty2, delay2, images2, duration2)

Builds the camera command string based on the parameters.
"""
function create_camera_command(freq, ch1Enabled, duty1, delay1, images1, duration1,
                               ch2Enabled, duty2, delay2, images2, duration2)
    cmd = ""
    if ch1Enabled || ch2Enabled
        cmd *= "S 3 S $images1\n"
        cmd *= "S 4 S $images2\n"

        cmd *= "A 3 WT $(round(Int, (100 - duration1) * duty1/images1/freq * 50))\n"
        cmd *= "A 3 DP 1 $freq $(round(Int, duration1 * duty1/images1/100))\n"
        cmd *= "A 3 WT $(round(Int, (100 + duration1) * duty1/images1/freq * 50))\n"

        cmd *= "A 4 WT $(round(Int, (100 - duration2) * duty2/images2/freq * 50))\n"
        cmd *= "A 4 DP 2 $freq $(round(Int, duration2 * duty2/images2/100))\n"
        cmd *= "A 4 WT $(round(Int, (100 + duration2) * duty2/images2/freq * 50))\n"

        # cmd_ch1 = ""
        # cmd_ch2 = ""
        # for _ in 1:images1
        #     cmd_ch1 *= "A 2 WT $(round(Int, (100 - duration1) * duty1/images1/freq * 50))\n"
        #     cmd_ch1 *= "A 2 DP 1 $freq $(duration1 * duty1/images1/100)\n"
        #     cmd_ch1 *= "A 2 WT $(round(Int, (100 + duration1) * duty1/images1/freq * 50))\n"
        # end
        # for _ in 1:images2
        #     cmd_ch2 *= "A 2 WT $(round(Int, (100 - duration2) * duty2/images2/freq * 50))\n"
        #     cmd_ch2 *= "A 2 DP 2 $freq $(duration2 * duty2/images2/100)\n"
        #     cmd_ch2 *= "A 2 WT $(round(Int, (100 + duration2) * duty2/images2/freq * 50))\n"
        # end
    
        if delay1 < delay2
            cmd *= "A 2 WT $(round(Int, delay1/freq * 10000))\n"
            if ch1Enabled
                cmd *= "A 2 TG 3\n"
            end
            cmd *= "A 2 WT $(round(Int, (delay2 - delay1)/freq * 10000))\n"
            if ch2Enabled
                cmd *= "A 2 TG 4\n"
            end
            cmd *= "A 2 WT $(round(Int, (100 - delay2)/freq * 10000))\n"
        else
            cmd *= "A 2 WT $(round(Int, delay2/freq * 10000))\n"
            if ch2Enabled
                cmd *= "A 2 TG 4\n"
            end
            cmd *= "A 2 WT $(round(Int, (delay1 - delay2)/freq * 10000))\n"
            if ch1Enabled
                cmd *= "A 2 TG 3\n"
            end
            cmd *= "A 2 WT $(round(Int, (100 - delay1)/freq * 10000))\n"
        end
    else
        println("No command available")
    end
    return cmd
end

"""
    create_command(freq, ch1Enabled, duty1, delay1, intensity1, images1, duration1,
                   ch2Enabled, duty2, delay2, intensity2, images2, duration2)

Builds the full command by concatenating system, LED (currently empty), and camera commands.
"""
function create_command(freq,
                        ch1Enabled, duty1, delay1, intensity1, images1, duration1,
                        ch2Enabled, duty2, delay2, intensity2, images2, duration2)
    # LED command is currently disabled.
    # LED_cmd = ""
    LED_cmd = create_LED_command(freq, ch1Enabled, duty1, delay1, intensity1, ch2Enabled, duty2, delay2, intensity2)
    # camera_cmd = ""
    camera_cmd = create_camera_command(freq, ch1Enabled, duty1, delay1, images1, duration1, ch2Enabled, duty2, delay2, images2, duration2)

    fullCmd = "C Z\nC R\n"

    if !isempty(LED_cmd)
        fullCmd *= "S 1 S 1\n" * LED_cmd
    end

    if !isempty(camera_cmd)
        fullCmd *= "S 2 S 1\n" * camera_cmd
    end

    if !isempty(LED_cmd)
        fullCmd *= "A 0 TG 1\n"
    end

    if !isempty(camera_cmd)
        fullCmd *= "A 0 TG 2\n"
    end

    # fullCmd = "C Z\nC R\nS 2 S 1\n" * LED_cmd * camera_cmd * "A 0 TG 2\n"
    # fullCmd = "C Z\nC R\nS 1 S 1\nS 2 S 1\n" * LED_cmd * camera_cmd * "A 0 TG 1\nA 0 TG 2\n"
    return fullCmd
end

########################################################################
# GUI & PLOTTING SETUP
########################################################################

# Set global theme
set_theme!(backgroundcolor = :gray95)

# Create main figure window
fig = Figure(size = (700, 800))

# Create layout grids
fig[1, 1] = plotGrid = GridLayout(tellwidth = true)
fig[2, 1] = controlGrid = GridLayout(tellwidth = false)
controlGrid[1, 1] = buttonGrid = GridLayout(tellwidth = false)
controlGrid[1, 2] = sliderGrid = GridLayout(tellwidth = false)
buttonGrid[2, 1] = freqGrid = GridLayout(tellwidth = false)
buttonGrid[2, 2] = fpsGrid  = GridLayout(tellwidth = false)

# Set column sizes and gaps
colsize!(controlGrid, 1, 200)
colgap!(controlGrid, 1, 64)

# Set colors
dark_gray = RGBf(0.4, 0.4, 0.4)
mid_gray = RGBf(0.6, 0.6, 0.6)
light_gray = RGBf(0.8, 0.8, 0.8)
off_white = RGBf(0.9, 0.9, 0.9)

# -------------------------------
# Frequency Controls
# -------------------------------

frequency_label = Label(buttonGrid[1, 1],
    fontsize = 16,
    justification = :right,
    text = "Frequency [Hz] :")

frequency_textbox = Textbox(freqGrid[1:2, 1],
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
frequency_up = Button(freqGrid[1, 2],
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

frequency_down = Button(freqGrid[2, 2],
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
Box(freqGrid[1:2, 1:2],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)

# -------------------------------
# FPS Controls
# -------------------------------

fps_label = Label(buttonGrid[1, 2],
    fontsize = 16,
    justification = :right,
    text = "Max fps:")

fps_textbox = Textbox(fpsGrid[1:2, 1],
    borderwidth = 0,
    boxcolor = :white,
    boxcolor_focused = :white,
    boxcolor_focused_invalid = RGBf(1, 0.9, 0.9),
    boxcolor_hover = off_white,
    cornerradius = 0,
    height = 32,
    placeholder = " ",
    validator = UInt16,
    width = 96)
fps_textbox.displayed_string[] = "1000"
fps_textbox.stored_string[] = "1000"

# FPS Up/Down Buttons
fps_up = Button(fpsGrid[1, 2],
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

fps_down = Button(fpsGrid[2, 2],
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
Box(fpsGrid[1:2, 1:2],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)

# Set gap adjustments
colgap!(freqGrid, 1, Fixed(-1))
rowgap!(freqGrid, 1, Fixed(0))
colgap!(fpsGrid, 1, Fixed(-1))
rowgap!(fpsGrid, 1, Fixed(0))
colgap!(buttonGrid, 1, Fixed(18))

# -------------------------------
# Start / Stop and Port Controls
# -------------------------------

start = Button(buttonGrid[3, 1],
    buttoncolor = light_gray,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = mid_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 16,
    height = 32,
    label = "START",
    labelcolor_active = :black,
    strokecolor = :black,
    strokewidth = 2,
    valign = :bottom,
    width = 128)

stop = Button(buttonGrid[3, 2],
    buttoncolor = light_gray,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = mid_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 16,
    height = 32,
    label = "STOP",
    labelcolor_active = :black,
    strokecolor = :black,
    strokewidth = 2,
    valign = :bottom,
    width = 128)

portname_textbox = Textbox(buttonGrid[4, 1],
    bordercolor = :black,
    bordercolor_focused = :black,
    bordercolor_hover = :black,
    borderwidth = 1,
    boxcolor = :white,
    boxcolor_focused = :white,
    boxcolor_focused_invalid = RGBf(1, 0.9, 0.9),
    boxcolor_hover = off_white,
    cornerradius = 0,
    height = 32,
    placeholder = "Port name...",
    width = 128)
portname_textbox.displayed_string[] = "usbmodem102"
portname_textbox.stored_string[] = "usbmodem102"

connect = Button(buttonGrid[4, 2],
    buttoncolor = light_gray,
    buttoncolor_active = dark_gray,
    buttoncolor_hover = mid_gray,
    cornerradius = 0,
    font = :bold,
    fontsize = 16,
    height = 32,
    label = "Connect",
    labelcolor_active = :black,
    strokecolor = :black,
    strokewidth = 2,
    valign = :bottom,
    width = 128)

monitor_box = Box(buttonGrid[5, 1:2],
    color = :white,
    cornerradius = 0,
    height = 296,
    width = 274)
global monitor = Label(buttonGrid[5, 1:2], "\n",
    height = 120,
    justification = :left,
    padding = (8, 8, 8, 8),
    tellheight = false,
    tellwidth = false,
    valign = :bottom,
    width = 274)

########################################################################
# SLIDER CONTROLS FOR CHANNELS
########################################################################

# Channel 1 Controls
ch1_check = Toggle(sliderGrid[1, 1][1, 1],
    active = true,
    buttoncolor = dark_gray,
    framecolor_active = mid_gray,
    framecolor_inactive = light_gray,
    rimfraction = 0.24)

ch1_label = Label(sliderGrid[1, 1][1, 2],
    font = :bold,
    fontsize = 18,
    text = "Channel 1")

ch1_menu = Menu(sliderGrid[1, 1][1, 3],
    cell_color_active = light_gray,
    cell_color_hover = off_white,
    cell_color_inactive_even = :white,
    cell_color_inactive_odd = :white,
    default = "Square",
    dropdown_arrow_color = :black,
    height = 34,
    options = ["Square", "Sawtooth", "R-sawtooth", "Triangle", "R-triangle", "Exponential"],
    selection_cell_color_inactive = :white,
    width = 128)

Box(sliderGrid[1, 1][1, 3],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)

ch1_sliders = SliderGrid(sliderGrid[1, 1][2, :],
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Duty cycle",
     linewidth = 14,
     range = 1:1:100,
     startvalue = 40,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Delay",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 5,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Intensity",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 100,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = " {:d} (10000 fps)",
     label = "Images",
     linewidth = 14,
     range = 0:1:40,
     startvalue = 5,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} % (999 ms)",
     label = "Pulse\nduration",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 50,
     tellheight = false),
    width = 350)
ch1_sliders.layout.content[12].content.text = "5 (125 fps)"
ch1_sliders.layout.content[15].content.text = "50 % (4.0 ms)"

ch1_visible = Box(sliderGrid[1, 1][2, :],
    color = RGBAf(0.95, 0.95, 0.95, 0.5),
    cornerradius = 0,
    height = 180,
    strokevisible = false,
    visible = false,
    width = 360)

# Channel 2 Controls
ch2_check = Toggle(sliderGrid[3, 1][1, 1],
    active = true,
    buttoncolor = dark_gray,
    framecolor_active = mid_gray,
    framecolor_inactive = light_gray,
    rimfraction = 0.24)

ch2_label = Label(sliderGrid[3, 1][1, 2],
    font = :bold,
    fontsize = 18,
    text = "Channel 2")

ch2_menu = Menu(sliderGrid[3, 1][1, 3],
    cell_color_active = light_gray,
    cell_color_hover = off_white,
    cell_color_inactive_even = :white,
    cell_color_inactive_odd = :white,
    default = "Square",
    dropdown_arrow_color = :black,
    height = 34,
    options = ["Square", "Sawtooth", "R-sawtooth", "Triangle", "R-triangle", "Exponential"],
    selection_cell_color_inactive = :white,
    width = 128)

Box(sliderGrid[3, 1][1, 3],
    color = :transparent,
    cornerradius = 0,
    height = 34,
    width = 128)

ch2_sliders = SliderGrid(sliderGrid[3, 1][2, :],
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Duty cycle",
     linewidth = 14,
     range = 1:1:100,
     startvalue = 40,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Delay",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 55,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} %",
     label = "Intensity",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 100,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = " {:d} (10000 fps)",
     label = "Images",
     linewidth = 14,
     range = 0:1:40,
     startvalue = 2,
     tellheight = false),
    (color_active = dark_gray,
     color_active_dimmed = mid_gray,
     color_inactive = light_gray,
     format = "{:d} % (999 ms)",
     label = "Pulse\nduration",
     linewidth = 14,
     range = 0:1:100,
     startvalue = 50,
     tellheight = false),
    width = 350)
ch2_sliders.layout.content[12].content.text = "2 (50 fps)"
ch2_sliders.layout.content[15].content.text = "50 % (10.0 ms)"

ch2_visible = Box(sliderGrid[3, 1][2, :],
    color = RGBAf(0.95, 0.95, 0.95, 0.5),
    cornerradius = 0,
    height = 180,
    strokevisible = false,
    visible = false,
    width = 360)

# Set fixed row gap for the slider grid
rowsize!(sliderGrid, 2, Fixed(20))

########################################################################
# SIGNAL & PLOT INITIALIZATION
########################################################################

numPeriods = 2
resolution = 20000
t_vals = range(0, numPeriods, resolution)
baseline = zeros(resolution)

global ch1_vector_gen = square_vector
global ch2_vector_gen = square_vector

global LED_1 = Observable(zeros(resolution))
LED_1[] = ch1_vector_gen(numPeriods, resolution, 40, 5, 100)

global LED_2 = Observable(zeros(resolution))
LED_2[] = ch2_vector_gen(numPeriods, resolution, 40, 55, 100)

global camera_1 = Observable(zeros(resolution))
camera_1[] = camera_vector(numPeriods, resolution, 40, 5, 5, 50) .- 10

global camera_2 = Observable(zeros(resolution))
camera_2[] = camera_vector(numPeriods, resolution, 40, 55, 2, 50) .- 10

global error = Observable(zeros(resolution))
error[] = camera_1[].*camera_2[]./100 .- 10

global serial = nothing
global connected = Observable(false)

global axis = plotGrid[1, 1] = Axis(fig,
    title = "Signal",
    ylabel = "Intensity [%]",
    titlesize = 16,
    aspect = 4,
    limits = (-0.1, numPeriods + 0.1, -9, 109),
    xticks = 0:1:numPeriods,
    xtickformat = "{:d} T")

# Plot bands and lines for signals
band!(axis, t_vals, LED_1, baseline, color=Makie.wong_colors()[1], alpha=0.25)
band!(axis, t_vals, LED_2, baseline, color=Makie.wong_colors()[2], alpha=0.25)
camera_1_band = band!(axis, t_vals, camera_1, baseline .- 10, color=:gray25, alpha=0.5)
camera_2_band = band!(axis, t_vals, camera_2, baseline .- 10, color=:gray25, alpha=0.5)
lines!(axis, t_vals, LED_1, color=Makie.wong_colors()[1])
lines!(axis, t_vals, LED_2, color=Makie.wong_colors()[2])
error_band = band!(axis, t_vals, error, baseline .- 10, color=:red2, alpha=0.5)
lines!(axis, t_vals, baseline, color=:black)

########################################################################
# EVENT HANDLERS & REACTIVE UPDATES
########################################################################

# Frequency textbox validation
on(frequency_textbox.stored_string) do str_val
    if parse(Float64, str_val) < 0.04
        frequency_textbox.displayed_string[] = "0.04"
        frequency_textbox.stored_string[] = "0.04"
    elseif parse(Float64, str_val) > 1_000_000
        frequency_textbox.displayed_string[] = "1000000.0"
        frequency_textbox.stored_string[] = "1000000.0"
    end
    max_images_ch1 = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, str_val)*ch1_sliders.sliders[1].value[]/100)
    max_images_ch2 = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, str_val)*ch2_sliders.sliders[1].value[]/100)
    ch1_sliders.sliders[4].range = 0:1:max_images_ch1
    ch2_sliders.sliders[4].range = 0:1:max_images_ch2
end

# Frequency up button
on(frequency_up.clicks) do _
    orig_freq = parse(Float64, frequency_textbox.stored_string[])
    freq_new = ceil(Int, orig_freq/(10^floor(Int, log10(orig_freq))) + 1e-9) * 10^floor(Int, log10(orig_freq))
    if freq_new < 0.04
        freq_new = 0.04
    elseif freq_new > 1_000_000
        freq_new = 1_000_000
    end
    frequency_textbox.displayed_string[] = string(freq_new)*".0"
    frequency_textbox.stored_string[] = string(freq_new)*".0"
end

# Frequency down button
on(frequency_down.clicks) do _
    orig_freq = parse(Float64, frequency_textbox.stored_string[])
    if log10(orig_freq) % 1 == 0
        freq_new = floor(Int, orig_freq/(10^floor(Int, log10(orig_freq)-1)) - 1e-9) * 10^floor(Int, log10(orig_freq)-1)
    else
        freq_new = floor(Int, orig_freq/(10^floor(Int, log10(orig_freq))) - 1e-9) * 10^floor(Int, log10(orig_freq))
    end
    if freq_new < 0.04
        freq_new = 0.04
    elseif freq_new > 1_000_000
        freq_new = 1_000_000
    end
    frequency_textbox.displayed_string[] = string(freq_new)*".0"
    frequency_textbox.stored_string[] = string(freq_new)*".0"
end

# FPS textbox validation
on(fps_textbox.stored_string) do str_val
    max_images_ch1 = floor(Int, parse(Int, str_val)/parse(Float64, frequency_textbox.stored_string[])*ch1_sliders.sliders[1].value[]/100)
    max_images_ch2 = floor(Int, parse(Int, str_val)/parse(Float64, frequency_textbox.stored_string[])*ch2_sliders.sliders[1].value[]/100)
    ch1_sliders.sliders[4].range = 0:1:max_images_ch1
    ch2_sliders.sliders[4].range = 0:1:max_images_ch2
end

# FPS up button
on(fps_up.clicks) do _
    orig_fps = parse(Float64, fps_textbox.stored_string[])
    fps_new = ceil(Int, orig_fps/(10^floor(Int, log10(orig_fps))) + 1e-9) * 10^floor(Int, log10(orig_fps))
    fps_textbox.displayed_string[] = string(fps_new)
    fps_textbox.stored_string[] = string(fps_new)
end

# FPS down button
on(fps_down.clicks) do _
    orig_fps = parse(Float64, fps_textbox.stored_string[])
    if log10(orig_fps) % 1 == 0
        fps_new = floor(Int, orig_fps/(10^floor(Int, log10(orig_fps)-1)) - 1e-9) * 10^floor(Int, log10(orig_fps)-1)
    else
        fps_new = floor(Int, orig_fps/(10^floor(Int, log10(orig_fps))) - 1e-9) * 10^floor(Int, log10(orig_fps))
    end
    fps_textbox.displayed_string[] = string(fps_new)
    fps_textbox.stored_string[] = string(fps_new)
end

# Start button – build command and send via serial
on(start.clicks) do _
    cmd = create_command(round(parse(Float64, frequency_textbox.displayed_string[]), digits=2),
            ch1_check.active[], ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[],
            ch1_sliders.sliders[3].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[],
            ch2_check.active[], ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[],
            ch2_sliders.sliders[3].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
    println(cmd)
    if connected[]
        if error[] < zeros(resolution)
            # Rebuild command and send
            cmd = create_command(round(parse(Float64, frequency_textbox.displayed_string[]), digits=2),
                ch1_check.active[], ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[],
                ch1_sliders.sliders[3].value[], ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[],
                ch2_check.active[], ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[],
                ch2_sliders.sliders[3].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[])
            send_command(serial, cmd)
        else
            add_to_monitor("Camera superposition error")
            println("Camera superposition error")
        end
    else
        add_to_monitor("Not connected")
        println("Not connected")
    end
end

# Stop button – send stop command
on(stop.clicks) do _
    if connected[]
        send_command(serial, "C Z\n")
    else
        add_to_monitor("Not connected")
        println("Not connected")
    end
end

# Connect/Disconnect button
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

########################################################################
# SLIDER EVENT HANDLERS FOR CHANNEL UPDATES
########################################################################

# Channel 1 – Duty cycle slider handler
on(ch1_sliders.sliders[1].value) do val
    if ch1_check.active[]
        LED_1[] = ch1_vector_gen(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                             ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[])
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                    ch1_sliders.sliders[2].value[], 120) .- 10
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[],
                                       ch1_sliders.sliders[5].value[]) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end

    max_images = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, frequency_textbox.stored_string[])*val/100)
    ch1_sliders.sliders[4].range = 0:1:max_images

    # Update fps text display for channel 1
    fps_calc = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/
                     parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
    if fps_calc >= 1
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps_calc ms)"
    else
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
    end
end

# Channel 1 – Delay slider handler
on(ch1_sliders.sliders[2].value) do val
    if ch1_check.active[]
        LED_1[] = ch1_vector_gen(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                             val, ch1_sliders.sliders[3].value[])
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                    val, 120) .- 10
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       val, ch1_sliders.sliders[4].value[], ch1_sliders.sliders[5].value[]) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
end

# Channel 1 – Intensity slider handler
on(ch1_sliders.sliders[3].value) do val
    if ch1_check.active[]
        LED_1[] = ch1_vector_gen(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                             ch1_sliders.sliders[2].value[], val)
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
end

# Channel 1 – Images slider handler
on(ch1_sliders.sliders[4].value) do val
    if ch1_check.active[]
        if val/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                    ch1_sliders.sliders[2].value[], 120) .- 10
            camera_1_band.alpha = ch1_sliders.sliders[5].value[]/200
            error_band.alpha = ch1_sliders.sliders[5].value[]/200
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       ch1_sliders.sliders[2].value[], val, ch1_sliders.sliders[5].value[]) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        ch1_sliders.layout.content[12].content.text = "$val ($(round(Int, 100*val*parse(Float64, frequency_textbox.stored_string[])/ch1_sliders.sliders[1].value[])) fps)"
        fps_calc = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/
                         parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
        if fps_calc >= 1
            ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps_calc ms)"
        else
            ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
        end
        
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
end

# Channel 1 – Pulse duration slider handler
on(ch1_sliders.sliders[5].value) do val
    if ch1_check.active[]
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                    ch1_sliders.sliders[2].value[], 120) .- 10
            camera_1_band.alpha = val/200
            error_band.alpha = val/200
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[], val) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
    fps_calc = round(ch1_sliders.sliders[1].value[]*ch1_sliders.sliders[5].value[]/ch1_sliders.sliders[4].value[]/
                     parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
    if fps_calc >= 1
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($fps_calc ms)"
    else
        ch1_sliders.layout.content[15].content.text = "$(ch1_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
    end
end

# Channel 2 – Duty cycle slider handler
on(ch2_sliders.sliders[1].value) do val
    if ch2_check.active[]
        LED_2[] = ch2_vector_gen(numPeriods, resolution, val, ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[])
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 120) .- 10
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(numPeriods, resolution, val, ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[]) .- 10
            camera_2_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end

    max_images = floor(Int, parse(Int, fps_textbox.stored_string[])/parse(Float64, frequency_textbox.stored_string[])*val/100)
    ch2_sliders.sliders[4].range = 0:1:max_images

    fps_calc = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/
                     parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
    if fps_calc >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps_calc ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
    end
end

# Channel 2 – Delay slider handler
on(ch2_sliders.sliders[2].value) do val
    if ch2_check.active[]
        LED_2[] = ch2_vector_gen(numPeriods, resolution, ch2_sliders.sliders[1].value[], val, ch2_sliders.sliders[3].value[])
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], val, 120) .- 10
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], val, ch2_sliders.sliders[4].value[], ch2_sliders.sliders[5].value[]) .- 10
            camera_2_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
end

# Channel 2 – Intensity slider handler
on(ch2_sliders.sliders[3].value) do val
    if ch2_check.active[]
        LED_2[] = ch2_vector_gen(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], val)
        error[] = camera_1[].*camera_2[]./100 .- 10
    end
end

# Channel 2 – Images slider handler
on(ch2_sliders.sliders[4].value) do val
    if ch2_check.active[]
        if val/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 120) .- 10
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], val, ch2_sliders.sliders[5].value[]) .- 10
            camera_2_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end

    ch2_sliders.layout.content[12].content.text = "$val ($(round(Int, 100*val*parse(Float64, frequency_textbox.stored_string[])/ch2_sliders.sliders[1].value[])) fps)"
    fps_calc = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/
                     parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
    if fps_calc >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps_calc ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
    end
end

# Channel 2 – Pulse duration slider handler
on(ch2_sliders.sliders[5].value) do val
    if ch2_check.active[]
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 120) .- 10
            camera_2_band.alpha = ch2_sliders.sliders[5].value[]/200
            error_band.alpha = ch2_sliders.sliders[5].value[]/200
        else
            camera_2[] = camera_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[], val) .- 10
            camera_2_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10
    end

    fps_calc = round(ch2_sliders.sliders[1].value[]*ch2_sliders.sliders[5].value[]/ch2_sliders.sliders[4].value[]/
                     parse(Float64, frequency_textbox.stored_string[])/10, sigdigits=2)
    if fps_calc >= 1
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($fps_calc ms)"
    else
        ch2_sliders.layout.content[15].content.text = "$(ch2_sliders.sliders[5].value[]) % ($(round(Int, 1000*fps_calc)) μs)"
    end
end

# Channel 1 – Enable/Disable checkbox handler
on(ch1_check.active) do active
    if active
        LED_1[] = ch1_vector_gen(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                             ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[])
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], 120) .- 10
            camera_1_band.alpha = 0.25
            error_band.alpha = 0.25
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[],
                                       ch1_sliders.sliders[5].value[]) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch1_visible.visible[] = false
    else
        LED_1[] = zeros(resolution)
        camera_1[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch1_visible.visible[] = true
    end
end

# Channel 2 – Enable/Disable checkbox handler
on(ch2_check.active) do active
    if active
        LED_2[] = ch2_vector_gen(numPeriods, resolution, ch2_sliders.sliders[1].value[],
                             ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[])
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 120) .- 10
            camera_2_band.alpha = 0.25
        else
            camera_2[] = camera_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[],
                                       ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[],
                                       ch2_sliders.sliders[5].value[]) .- 10
            camera_2_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch2_visible.visible[] = false
    else
        LED_2[] = zeros(resolution)
        camera_2[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch2_visible.visible[] = true
    end
end

on(ch1_menu.selection) do type
    if type == "Square"
        global ch1_vector_gen = square_vector
    elseif type == "Sawtooth"
        global ch1_vector_gen = saw_vector
    elseif type == "R-sawtooth"
        global ch1_vector_gen = rsaw_vector
    elseif type == "Triangle"
        global ch1_vector_gen = tri_vector
    elseif type == "R-triangle"
        global ch1_vector_gen = rtri_vector
    elseif type == "Exponential"
        global ch1_vector_gen = exp_vector
    else
        println("Vector type not implemented")
    end

    if ch1_check.active[]
        LED_1[] = ch1_vector_gen(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                             ch1_sliders.sliders[2].value[], ch1_sliders.sliders[3].value[])
        if ch1_sliders.sliders[4].value[]/ch1_sliders.sliders[1].value[] >= 2/numPeriods
            camera_1[] = square_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[], ch1_sliders.sliders[2].value[], 120) .- 10
            camera_1_band.alpha = 0.25
            error_band.alpha = 0.25
        else
            camera_1[] = camera_vector(numPeriods, resolution, ch1_sliders.sliders[1].value[],
                                       ch1_sliders.sliders[2].value[], ch1_sliders.sliders[4].value[],
                                       ch1_sliders.sliders[5].value[]) .- 10
            camera_1_band.alpha = 0.5
            error_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch1_visible.visible[] = false
    else
        LED_1[] = zeros(resolution)
        camera_1[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch1_visible.visible[] = true
    end
end

on(ch2_menu.selection) do type
    if type == "Square"
        global ch2_vector_gen = square_vector
    elseif type == "Sawtooth"
        global ch2_vector_gen = saw_vector
    elseif type == "R-sawtooth"
        global ch2_vector_gen = rsaw_vector
    elseif type == "Triangle"
        global ch2_vector_gen = tri_vector
    elseif type == "R-triangle"
        global ch2_vector_gen = rtri_vector
    elseif type == "Exponential"
        global ch2_vector_gen = exp_vector
    else
        println("Vector type not implemented")
    end

    if ch2_check.active[]
        LED_2[] = ch2_vector_gen(numPeriods, resolution, ch2_sliders.sliders[1].value[],
                             ch2_sliders.sliders[2].value[], ch2_sliders.sliders[3].value[])
        if ch2_sliders.sliders[4].value[]/ch2_sliders.sliders[1].value[] >= 2/numPeriods
            camera_2[] = square_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[], ch2_sliders.sliders[2].value[], 120) .- 10
            camera_2_band.alpha = 0.25
        else
            camera_2[] = camera_vector(numPeriods, resolution, ch2_sliders.sliders[1].value[],
                                       ch2_sliders.sliders[2].value[], ch2_sliders.sliders[4].value[],
                                       ch2_sliders.sliders[5].value[]) .- 10
            camera_2_band.alpha = 0.5
        end
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch2_visible.visible[] = false
    else
        LED_2[] = zeros(resolution)
        camera_2[] = zeros(resolution)
        error[] = camera_1[].*camera_2[]./100 .- 10

        ch2_visible.visible[] = true
    end
end

########################################################################
# DISPLAY THE FIGURE
########################################################################
fig
