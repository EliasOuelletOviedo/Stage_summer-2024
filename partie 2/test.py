#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk, scrolledtext
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import serial
import time, math

###############################################################################
# SIGNAL GENERATION FUNCTIONS
###############################################################################
def pwm_vector(periods, resolution, duty_cycle, delay, intensity):
    """
    Generate a PWM vector for LED signals.
    """
    vector = np.zeros(resolution)
    num_ones = int(round(resolution * duty_cycle / 100 / periods))
    start_idx = int(round(resolution * delay / 100 / periods))
    for period in range(1, periods + 1):
        for i in range(num_ones):
            # Note: Python indexing is 0-based.
            idx = int(round((start_idx + i + (period + 1)/periods * resolution))) % resolution
            vector[idx] = intensity
    return vector

def camera_vector(periods, resolution, duty_cycle, delay, images, duration):
    """
    Generate a vector for camera trigger signals.
    """
    vector = np.zeros(resolution)
    num_ones  = int(round(duration * resolution * duty_cycle / 5000 / periods))
    num_zeros = int(round((100 - duration) * resolution * duty_cycle / 5000 / periods))
    start_idx = int(round(resolution * delay / 100 / periods))
    if images == 0:
        return vector
    pulse_duration = num_ones // images
    pulse_interval = num_zeros // images
    for period in range(1, periods + 1):
        for i in range(images):
            for j in range(pulse_duration // 2):
                idx = int(round((start_idx +
                                  pulse_interval/4 +
                                  i * (pulse_duration + pulse_interval)/2 +
                                  j +
                                  (period + 1)/periods * resolution))) % resolution
                vector[idx] = 100
    return vector

###############################################################################
# SERIAL COMMUNICATION FUNCTIONS
###############################################################################
def add_to_monitor(monitor_widget, msg):
    """Append a message to the monitor text widget (keeps last 17 lines)."""
    current = monitor_widget.get("1.0", tk.END).splitlines()
    current.append(msg.splitlines()[0])
    MAX_LINES = 17
    if len(current) > MAX_LINES:
        current = current[-MAX_LINES:]
    monitor_widget.delete("1.0", tk.END)
    monitor_widget.insert(tk.END, "\n".join(current))

def connect_to_port(portname):
    """Attempts to connect via serial to '/dev/tty.<portname>'."""
    try:
        ser = serial.Serial(f"/dev/tty.{portname}", 115200, timeout=3)
        time.sleep(0.1)
        print("Connection successful")
        return ser
    except Exception as e:
        print("Connection error:", e)
        return None

def send_command(ser, command, monitor_widget):
    """Sends a command over serial and updates the monitor with the response."""
    print("Sending command:")
    print(command)
    try:
        ser.write(command.encode())
        time.sleep(0.1)
        response = ser.read_all().decode().strip()
    except Exception as e:
        response = f"Error during send: {e}"
    print(response)
    add_to_monitor(monitor_widget, response)
    return response

###############################################################################
# COMMAND CREATION FUNCTIONS
###############################################################################
def create_LED_command(frequency, ch1, duty_cycle_1, delay_1, intensity_1,
                       ch2, duty_cycle_2, delay_2, intensity_2):
    command = ""
    if ch1 or ch2:
        channel_1_cmd = f"A 1 AP 1 {frequency} {duty_cycle_1} 0 {int(intensity_1 * 100)}\n"
        channel_2_cmd = f"A 1 AP 2 {frequency} {duty_cycle_2} 0 {int(intensity_2 * 100)}\n"
        if delay_1 < delay_2:
            command += f"A 1 WT {int(round(delay_1/frequency * 10000))}\n"
            if ch1:
                command += channel_1_cmd
            command += f"A 1 WT {int(round((delay_2 - delay_1)/frequency * 10000))}\n"
            if ch2:
                command += channel_2_cmd
            command += f"A 1 WT {int(round((100 - delay_2)/frequency * 10000))}\n"
        else:
            command += f"A 1 WT {int(round(delay_1/frequency * 10000))}\n"
            if ch2:
                command += channel_2_cmd
            command += f"A 1 WT {int(round((delay_1 - delay_2)/frequency * 10000))}\n"
            if ch1:
                command += channel_1_cmd
            command += f"A 1 WT {int(round((100 - delay_1)/frequency * 10000))}\n"
    else:
        print("No LED command available (both channels disabled)")
    return command

def create_camera_command(frequency, ch1, duty_cycle_1, delay_1, images_1, duration_1,
                          ch2, duty_cycle_2, delay_2, images_2, duration_2):
    command = ""
    if ch1 or ch2:
        channel_1_cmd = ""
        channel_2_cmd = ""
        for _ in range(images_1):
            channel_1_cmd += f"A 1 WT {int(round((100 - duration_1)*duty_cycle_1/images_1/frequency*50))}\n"
            channel_1_cmd += f"A 1 DP 1 {frequency} {duration_1*duty_cycle_1/images_1/100}\n"
            channel_1_cmd += f"A 1 WT {int(round((100 + duration_1)*duty_cycle_1/images_1/frequency*50))}\n"
        for _ in range(images_2):
            channel_2_cmd += f"A 1 WT {int(round((100 - duration_2)*duty_cycle_2/images_2/frequency*50))}\n"
            channel_2_cmd += f"A 1 DP 2 {frequency} {duration_2*duty_cycle_2/images_2/100}\n"
            channel_2_cmd += f"A 1 WT {int(round((100 + duration_2)*duty_cycle_2/images_2/frequency*50))}\n"
        if delay_1 < delay_2:
            command += f"A 1 WT {int(round(delay_1/frequency * 10000))}\n"
            if ch1:
                command += channel_1_cmd
            command += f"A 1 WT {int(round((delay_2 - delay_1 - duty_cycle_1)/frequency * 10000))}\n"
            if ch2:
                command += channel_2_cmd
            command += f"A 1 WT {int(round((100 - delay_2 - duty_cycle_2)/frequency * 10000))}\n"
        else:
            command += f"A 1 WT {int(round(delay_1/frequency * 10000))}\n"
            if ch2:
                command += channel_2_cmd
            command += f"A 1 WT {int(round((delay_1 - delay_2 - duty_cycle_2)/frequency * 10000))}\n"
            if ch1:
                command += channel_1_cmd
            command += f"A 1 WT {int(round((100 - delay_1 - duty_cycle_1)/frequency * 10000))}\n"
    else:
        print("No camera command available (both channels disabled)")
    return command

def create_command(frequency,
                   ch1_enabled, duty_cycle_1, delay_1, intensity_1, images_1, duration_1,
                   ch2_enabled, duty_cycle_2, delay_2, intensity_2, images_2, duration_2):
    # The LED command is currently not active (empty) as in the original code.
    LED_command = ""
    camera_command = create_camera_command(frequency,
                                           ch1_enabled, duty_cycle_1, delay_1, images_1, duration_1,
                                           ch2_enabled, duty_cycle_2, delay_2, images_2, duration_2)
    command = "C Z\nC R\nS 1 S 1\n" + LED_command + camera_command + "A 0 TG 1\n"
    return command

###############################################################################
# TKINTER GUI APPLICATION
###############################################################################
class App:
    def __init__(self, root):
        self.root = root
        root.title("GUI Translator")
        
        # Global state variables
        self.periods = 2
        self.resolution = 20000
        self.serial_conn = None
        self.connected = False
        
        # Signal variables (will be NumPy arrays)
        self.LED_1 = pwm_vector(self.periods, self.resolution, 40, 5, 100)
        self.LED_2 = pwm_vector(self.periods, self.resolution, 40, 55, 100)
        self.camera_1 = camera_vector(self.periods, self.resolution, 40, 5, 5, 50)
        self.camera_2 = camera_vector(self.periods, self.resolution, 40, 55, 2, 50)
        self.error_signal = self.camera_1 * self.camera_2 / 100
        
        # Create frames for layout
        self.left_frame = ttk.Frame(root)
        self.left_frame.grid(row=0, column=0, sticky="nsew")
        self.right_frame = ttk.Frame(root)
        self.right_frame.grid(row=0, column=1, sticky="nsew")
        
        # Configure grid weights
        root.columnconfigure(0, weight=3)
        root.columnconfigure(1, weight=2)
        root.rowconfigure(0, weight=1)
        
        self.create_plot_area()
        self.create_controls()
    
    def create_plot_area(self):
        """Set up the Matplotlib plot inside a Tkinter canvas."""
        self.fig, self.ax = plt.subplots(figsize=(7, 4))
        self.ax.set_title("Signal")
        self.ax.set_ylabel("Intensity [%]")
        self.ax.set_xlim(-0.1, self.periods + 0.1)
        self.ax.set_ylim(-10, 110)
        self.ax.set_xticks(np.arange(0, self.periods+1, 1))
        self.line_LED1, = self.ax.plot(self.LED_1, color="C0", label="LED 1")
        self.line_LED2, = self.ax.plot(self.LED_2, color="C1", label="LED 2")
        self.line_camera1, = self.ax.plot(self.camera_1, color="black", alpha=0.5, label="Camera 1")
        self.line_camera2, = self.ax.plot(self.camera_2, color="black", alpha=0.5, label="Camera 2")
        self.ax.legend()
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.left_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def create_controls(self):
        """Create the right-side controls (labels, entries, sliders, buttons, monitor)."""
        # Top frame for frequency and fps controls.
        control_frame = ttk.Frame(self.right_frame)
        control_frame.pack(pady=5)
        
        # Frequency controls
        ttk.Label(control_frame, text="Frequency [Hz]:").grid(row=0, column=0, sticky="e")
        self.freq_var = tk.StringVar(value="10.0")
        self.freq_entry = ttk.Entry(control_frame, textvariable=self.freq_var, width=10)
        self.freq_entry.grid(row=0, column=1)
        self.freq_up = ttk.Button(control_frame, text="▲", command=self.freq_up_clicked)
        self.freq_up.grid(row=0, column=2)
        self.freq_down = ttk.Button(control_frame, text="▼", command=self.freq_down_clicked)
        self.freq_down.grid(row=0, column=3)
        
        # FPS controls
        ttk.Label(control_frame, text="Max fps:").grid(row=1, column=0, sticky="e")
        self.fps_var = tk.StringVar(value="1000")
        self.fps_entry = ttk.Entry(control_frame, textvariable=self.fps_var, width=10)
        self.fps_entry.grid(row=1, column=1)
        self.fps_up = ttk.Button(control_frame, text="▲", command=self.fps_up_clicked)
        self.fps_up.grid(row=1, column=2)
        self.fps_down = ttk.Button(control_frame, text="▼", command=self.fps_down_clicked)
        self.fps_down.grid(row=1, column=3)
        
        # Start/Stop and Connect buttons
        button_frame = ttk.Frame(self.right_frame)
        button_frame.pack(pady=5)
        self.start_btn = ttk.Button(button_frame, text="START", command=self.start_clicked)
        self.start_btn.grid(row=0, column=0, padx=5)
        self.stop_btn = ttk.Button(button_frame, text="STOP", command=self.stop_clicked)
        self.stop_btn.grid(row=0, column=1, padx=5)
        
        # Port name and Connect/Disconnect
        port_frame = ttk.Frame(self.right_frame)
        port_frame.pack(pady=5)
        ttk.Label(port_frame, text="Port name:").grid(row=0, column=0, sticky="e")
        self.port_var = tk.StringVar(value="usbmodem2303")
        self.port_entry = ttk.Entry(port_frame, textvariable=self.port_var, width=15)
        self.port_entry.grid(row=0, column=1)
        self.connect_btn = ttk.Button(port_frame, text="Connect", command=self.connect_clicked)
        self.connect_btn.grid(row=0, column=2, padx=5)
        
        # Monitor (scrolled text)
        self.monitor = scrolledtext.ScrolledText(self.right_frame, width=40, height=10)
        self.monitor.pack(pady=5)
        
        # Channel controls – For each channel, add a checkbutton and sliders.
        # Channel 1
        ch1_frame = ttk.LabelFrame(self.right_frame, text="Channel 1")
        ch1_frame.pack(fill="x", padx=5, pady=5)
        self.ch1_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(ch1_frame, text="Enabled", variable=self.ch1_enabled,
                        command=self.update_signals).grid(row=0, column=0, sticky="w")
        ttk.Label(ch1_frame, text="Duty cycle:").grid(row=1, column=0, sticky="e")
        self.ch1_duty = tk.Scale(ch1_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                 command=lambda v: self.update_signals())
        self.ch1_duty.set(40)
        self.ch1_duty.grid(row=1, column=1)
        ttk.Label(ch1_frame, text="Delay:").grid(row=2, column=0, sticky="e")
        self.ch1_delay = tk.Scale(ch1_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                  command=lambda v: self.update_signals())
        self.ch1_delay.set(5)
        self.ch1_delay.grid(row=2, column=1)
        ttk.Label(ch1_frame, text="Intensity:").grid(row=3, column=0, sticky="e")
        self.ch1_intensity = tk.Scale(ch1_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                      command=lambda v: self.update_signals())
        self.ch1_intensity.set(100)
        self.ch1_intensity.grid(row=3, column=1)
        ttk.Label(ch1_frame, text="Images:").grid(row=4, column=0, sticky="e")
        self.ch1_images = tk.Scale(ch1_frame, from_=0, to=40, orient=tk.HORIZONTAL,
                                   command=lambda v: self.update_signals())
        self.ch1_images.set(5)
        self.ch1_images.grid(row=4, column=1)
        ttk.Label(ch1_frame, text="Pulse duration:").grid(row=5, column=0, sticky="e")
        self.ch1_pulse = tk.Scale(ch1_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                  command=lambda v: self.update_signals())
        self.ch1_pulse.set(50)
        self.ch1_pulse.grid(row=5, column=1)
        
        # Channel 2
        ch2_frame = ttk.LabelFrame(self.right_frame, text="Channel 2")
        ch2_frame.pack(fill="x", padx=5, pady=5)
        self.ch2_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(ch2_frame, text="Enabled", variable=self.ch2_enabled,
                        command=self.update_signals).grid(row=0, column=0, sticky="w")
        ttk.Label(ch2_frame, text="Duty cycle:").grid(row=1, column=0, sticky="e")
        self.ch2_duty = tk.Scale(ch2_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                 command=lambda v: self.update_signals())
        self.ch2_duty.set(40)
        self.ch2_duty.grid(row=1, column=1)
        ttk.Label(ch2_frame, text="Delay:").grid(row=2, column=0, sticky="e")
        self.ch2_delay = tk.Scale(ch2_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                  command=lambda v: self.update_signals())
        self.ch2_delay.set(55)
        self.ch2_delay.grid(row=2, column=1)
        ttk.Label(ch2_frame, text="Intensity:").grid(row=3, column=0, sticky="e")
        self.ch2_intensity = tk.Scale(ch2_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                      command=lambda v: self.update_signals())
        self.ch2_intensity.set(100)
        self.ch2_intensity.grid(row=3, column=1)
        ttk.Label(ch2_frame, text="Images:").grid(row=4, column=0, sticky="e")
        self.ch2_images = tk.Scale(ch2_frame, from_=0, to=40, orient=tk.HORIZONTAL,
                                   command=lambda v: self.update_signals())
        self.ch2_images.set(2)
        self.ch2_images.grid(row=4, column=1)
        ttk.Label(ch2_frame, text="Pulse duration:").grid(row=5, column=0, sticky="e")
        self.ch2_pulse = tk.Scale(ch2_frame, from_=0, to=100, orient=tk.HORIZONTAL,
                                  command=lambda v: self.update_signals())
        self.ch2_pulse.set(50)
        self.ch2_pulse.grid(row=5, column=1)
    
    def update_plot(self):
        """Recalculate signals and update the plot lines."""
        freq = float(self.freq_var.get())
        # Recalculate signals based on channel settings.
        if self.ch1_enabled.get():
            self.LED_1 = pwm_vector(self.periods, self.resolution, self.ch1_duty.get(), self.ch1_delay.get(), self.ch1_intensity.get())
            # For camera_1, if images/duty are large, use pwm_vector instead of camera_vector.
            if self.ch1_images.get() / self.ch1_duty.get() >= 2/self.periods:
                self.camera_1 = pwm_vector(self.periods, self.resolution, self.ch1_duty.get(), self.ch1_delay.get(), 100)
            else:
                self.camera_1 = camera_vector(self.periods, self.resolution, self.ch1_duty.get(), self.ch1_delay.get(), self.ch1_images.get(), self.ch1_pulse.get())
        else:
            self.LED_1 = np.zeros(self.resolution)
            self.camera_1 = np.zeros(self.resolution)
            
        if self.ch2_enabled.get():
            self.LED_2 = pwm_vector(self.periods, self.resolution, self.ch2_duty.get(), self.ch2_delay.get(), self.ch2_intensity.get())
            if self.ch2_images.get() / self.ch2_duty.get() >= 2/self.periods:
                self.camera_2 = pwm_vector(self.periods, self.resolution, self.ch2_duty.get(), self.ch2_delay.get(), 100)
            else:
                self.camera_2 = camera_vector(self.periods, self.resolution, self.ch2_duty.get(), self.ch2_delay.get(), self.ch2_images.get(), self.ch2_pulse.get())
        else:
            self.LED_2 = np.zeros(self.resolution)
            self.camera_2 = np.zeros(self.resolution)
            
        self.error_signal = self.camera_1 * self.camera_2 / 100
        
        # Update plot lines (using full resolution arrays; you might want to downsample for performance)
        self.line_LED1.set_ydata(self.LED_1)
        self.line_LED2.set_ydata(self.LED_2)
        self.line_camera1.set_ydata(self.camera_1)
        self.line_camera2.set_ydata(self.camera_2)
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
    
    def update_signals(self):
        self.update_plot()
    
    def freq_up_clicked(self):
        try:
            orig = float(self.freq_var.get())
        except ValueError:
            orig = 10.0
        factor = 10 ** math.floor(math.log10(orig))
        new_freq = math.ceil(orig / factor + 1e-9) * factor
        # clamp value:
        if new_freq < 0.04:
            new_freq = 0.04
        elif new_freq > 1_000_000:
            new_freq = 1_000_000
        self.freq_var.set(f"{new_freq:.1f}")
        self.update_signals()
    
    def freq_down_clicked(self):
        try:
            orig = float(self.freq_var.get())
        except ValueError:
            orig = 10.0
        # Use a similar approach as freq_up_clicked but with floor.
        if math.log10(orig) % 1 == 0:
            factor = 10 ** math.floor(math.log10(orig) - 1)
        else:
            factor = 10 ** math.floor(math.log10(orig))
        new_freq = math.floor(orig / factor - 1e-9) * factor
        if new_freq < 0.04:
            new_freq = 0.04
        elif new_freq > 1_000_000:
            new_freq = 1_000_000
        self.freq_var.set(f"{new_freq:.1f}")
        self.update_signals()
    
    def fps_up_clicked(self):
        try:
            orig = float(self.fps_var.get())
        except ValueError:
            orig = 1000
        factor = 10 ** math.floor(math.log10(orig))
        new_fps = math.ceil(orig / factor + 1e-9) * factor
        self.fps_var.set(str(int(new_fps)))
    
    def fps_down_clicked(self):
        try:
            orig = float(self.fps_var.get())
        except ValueError:
            orig = 1000
        if math.log10(orig) % 1 == 0:
            factor = 10 ** math.floor(math.log10(orig) - 1)
        else:
            factor = 10 ** math.floor(math.log10(orig))
        new_fps = math.floor(orig / factor - 1e-9) * factor
        self.fps_var.set(str(int(new_fps)))
    
    def start_clicked(self):
        """Build the command from the current settings and send it if connected."""
        try:
            freq = round(float(self.freq_var.get()), 2)
        except ValueError:
            freq = 10.0
        command = create_command(freq,
                                 self.ch1_enabled.get(), self.ch1_duty.get(), self.ch1_delay.get(), self.ch1_intensity.get(),
                                 self.ch1_images.get(), self.ch1_pulse.get(),
                                 self.ch2_enabled.get(), self.ch2_duty.get(), self.ch2_delay.get(), self.ch2_intensity.get(),
                                 self.ch2_images.get(), self.ch2_pulse.get())
        print("Command to send:")
        print(command)
        if self.connected and self.serial_conn:
            # For simplicity, we do not check for camera superposition error here.
            send_command(self.serial_conn, command, self.monitor)
        else:
            add_to_monitor(self.monitor, "Not connected")
    
    def stop_clicked(self):
        if self.connected and self.serial_conn:
            send_command(self.serial_conn, "C Z\n", self.monitor)
        else:
            add_to_monitor(self.monitor, "Not connected")
    
    def connect_clicked(self):
        if not self.connected:
            portname = self.port_var.get()
            self.serial_conn = connect_to_port(portname)
            if self.serial_conn:
                self.connected = True
                self.connect_btn.config(text="Disconnect")
                add_to_monitor(self.monitor, "Connection successful")
            else:
                add_to_monitor(self.monitor, "Connection error")
        else:
            try:
                if self.serial_conn:
                    self.serial_conn.close()
                self.connected = False
                self.connect_btn.config(text="Connect")
                add_to_monitor(self.monitor, "Disconnection successful")
            except Exception as e:
                add_to_monitor(self.monitor, f"Error disconnecting: {e}")

###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
