import serial as ser
import time

serial = ser.Serial(
    port="/dev/tty.usbmodem2303",
    baudrate=115200,
    bytesize=ser.EIGHTBITS,
    parity=ser.PARITY_NONE,
    stopbits=ser.STOPBITS_ONE,
    timeout=3
)

def send_command(command):
    serial.write(command.encode())

    response = serial.readline().decode().strip()
    return response

trigger = 1

frequency = 10
duty_cycle = 25

command = f"C Z\n"
command += f"C R\n"
command += f"S {trigger} S 1\n"

command += f"A {trigger} DP 1 {frequency} {duty_cycle}\n"
command += f"A {trigger} WT {500000/frequency}\n"
command += f"A {trigger} DP 2 {frequency} {duty_cycle}\n"

command += f"A 0 TG {trigger}\n"

print(command)

try:
    response = send_command(command)
    print("Response:", response)
except Exception as e:
    print("Error:", e)
finally:
    serial.close()
