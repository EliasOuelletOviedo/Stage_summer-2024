using LibSerialPort

try
    serial = LibSerialPort.open("/dev/tty.usbmodem2303", 115200)
    LibSerialPort.set_read_timeout(serial, 3)
catch e
    println(e)
end

sleep(1)

function send_command(serial, command)
    write(serial, command)
    response = strip(read(serial, String))
    return response
end

trigger = 1
frequency = 10
duty_cycle = 25
LEDs = 1

command = "C Z\n"
command *= "C R\n"
command *= "S $trigger S 1\n"
command *= "A $trigger DP 1 $frequency $duty_cycle\n"
command *= "A $trigger WT $(Int(500000 / frequency))\n"
command *= "A $trigger DP 2 $frequency $duty_cycle\n"
command *= "A 0 TG $trigger\n"

println(command)

try
    response = send_command(serial, command)
    println("Response: ", response)
catch e
    println("Error: ", e)
finally
    LibSerialPort.close(serial)
end
