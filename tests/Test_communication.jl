using LibSerialPort

# try
serial = LibSerialPort.open("/dev/tty.usbmodem102", 115200)
# LibSerialPort.set_read_timeout(serial, 3)
# catch e
    # println(e)
# end

# sleep(1)

function send_command(serial, command)
    write(serial, command)
    sleep(0.1)
    response = strip(read(serial, String))
    return response
end

function get_response(serial)
    response = strip(read(serial, String))
    return response
end

# trigger = 1
# frequency = 10
# duty_cycle = 25
# LEDs = 1

# command = "C Z\n"
# command *= "C R\n"
# command *= "S $trigger S 1\n"
# command *= "A $trigger DP 1 $frequency $duty_cycle\n"
# command *= "A $trigger WT $(Int(500000 / frequency))\n"
# command *= "A $trigger DP 2 $frequency $duty_cycle\n"
# command *= "A 0 TG $trigger\n"

# command = "C Z\n"
# command *= "C R\n"
# command *= "S 1 D 4 R"
# command *= "A 1 AW 1 10 -1000 1000"
# command *= "A 1 RS S 10 -1000 1000 -1000 1000 10"

# command *= "A 0 DP 1 10 1"

# wait_s = 50000

# wait = 100000

# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "S 2 S 1\n"

# command *= "A 1 DO 1 1\n"

# command *= "A 2 WT $wait\n"
# command *= "A 2 DO 1 0\n"
# command *= "A 2 WT $wait_s\n"
# command *= "A 2 DO 1 1\n"

# for x in 0:2
#     x_V = 1000*x/2 - 500
#     global command *= "A 1 AO 1 $x_V\n"
#     for y in 0:2
#         y_V = 1000*y/2 - 500
#         global command *= "A 1 AO 2 $y_V\n"
#         global command *= "A 1 WT $(wait + wait_s)\n"
#     end
# end

# command *= "A 0 TG 1\n"
# command *= "A 0 TG 2\n"


# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "A 1 TO 100 100 5000 500 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 TO 100 100 5000 -500 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 TO 100 100 5000 -500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 TO 100 100 5000 500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 0 TG 1\n"


# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 TO 200 200 5000 500 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 TO 200 200 5000 -500 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 TO 200 200 5000 -500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 TO 200 200 5000 500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 0 TG 1\n"

# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 AO 1 500\n"
# command *= "A 1 AO 2 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 AO 1 -500\n"
# command *= "A 1 AO 2 500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 AO 1 -500\n"
# command *= "A 1 AO 2 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 AO 1 500\n"
# command *= "A 1 AO 2 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait_s\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 0 TG 1\n"

# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "A 1 DO 1 1\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 DO 1 0\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 TO 100 100 5000 -500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 1 TO 100 100 5000 500 -500\n"
# command *= "A 1 WT $wait\n"
# command *= "A 0 TG 1\n"

# command = "C Z\n"

# for (n, i) in enumerate(-10:10)
#     # println(n, i)
#     global command *= "A 0 AB 1 $n $(i*100)\n"
# end

# command *= "A 0 AR 1 20 21\n"

# first = true
# n_x = 5
# n_y = 5
# r = 250
# m = 5
# f = 10

# command = "C Z\n"

# for x_V in range(-1000, 1000, n_x)
#     for y_V in range(-1000, 1000, n_y)
#         println(x_V, " ",y_V)
#         if first
#             global command *= "A 0 TO $r $(floor(Int, r/m)) 0.001 $(round(Int, y_V)) $(round(Int, x_V)) $(round(1/(n_x*n_y), digits=2)) B\n"
#             global first = false
#         else
#             global command *= "A 0 TO $r $(floor(Int, r/m)) 0.001 $(round(Int, y_V)) $(round(Int, x_V)) $(round(1/(n_x*n_y), digits=2)) A\n"
#         end
#     end
# end

# command *= "A 0 AA $f\n"
# command *= "A 0 DP 2 $(f*n_x*n_y) 10\n"


# command = "C Z\n"
# command *= "A 0 TO $r $(r/m) $(t/pi/m) 0 0\n"

f = 16.67
n = 512

# command = "C Z\n"
# command *= "S 1 D 1 R 1\n"
# command *= "S 2 D 1 R 1\n"

# command *= "A 1 AW 1 $(f) -1000 1000\n"
# command *= "A 2 AW 2 $(f/n) -1000 1000\n"

# command *= "A 0 DP 1 $f 10\n"
# command *= "A 0 DP 2 $(f/n) $(10)\n"



command = "C Z\n"
command *= "S 1 D 1 R 1\n"
command *= "S 2 D 1 R 1\n"

# command *= "A 1 WT $(100)\n"
command *= "A 1 AW 1 $f 1000 -1000\n"

command *= "A 2 AW 2 $(f/n) -1000 1000\n"
command *= "A 2 TD $(n-1)\n"

response = send_command(serial, command)
println("Response: ", response)

# for (m, i) in enumerate(range(-1000, 1000, n))
#     # println(n, i)
#     global command = "A 0 AB 2 $m $(round(Int, i))\n"
#     response = send_command(serial, command)
#     println("Response: ", response)
# end

# command = "A 2 AR 2 $(f/n) $n\n"

# command *= "A 2 TD $(n-1)\n"

# command *= "A 0 DP 1 $f 5\n"
# response = send_command(serial, command)
# println("Response: ", response)



# command = "C Z\n"
# command *= "S 1 S\n"
# command *= "A 1 AO 3 5000\n"
# command *= "A 1 WT 50\n"
# command *= "A 1 AO 3 0\n"
# command *= "A 1 WT 59950\n"

# command *= "S 2 D 1 R 1\n"
# command *= "A 2 DO 1 1\n"
# command *= "A 2 WT 50\n"
# command *= "A 2 DO 1 0\n"
# command *= "A 2 WT 59950\n"

# command *= "A 0 TG 1\n"
# send_command(serial, command)

# r = 200
# m = 2
# t = 100000
# dt = 10000
# command = "C Z\n"
# command *= "S 1 S 0\n"
# command *= "S 2 S 1\n"

# command *= "A 2 WT $t\n"
# command *= "A 2 AO 1 0\n"
# command *= "A 2 AO 2 0\n"
# command *= "A 2 DO 1 0\n"
# command *= "A 2 WT $dt\n"
# command *= "A 2 DO 1 1\n"

# command *= "A 1 TO $r $(r/m) $(5*t/pi/m) 0 0\n"
# response = send_command(serial, command)

# command = "A 1 TG 2\n"
# command *= "A 1 WT $(t + dt)\n"

# # command *= "A 1 DO 1 1\n"


# command *= "A 0 TG 1\n"




# println(command)

# try
# response = send_command(serial, command)
# println("Response: ", response)
# catch e
#     println("Error: ", e)
# finally
#     LibSerialPort.close(serial)
# end


# for i in 1:10
#     for v in V
#         # global command = "A 1 AO 1 $v\n"
#         send_command(serial, "C 1\nA 2 AO 2 $v\n")
#         sleep(1/f)
#     end
#     # for v in reverse(V)
#     #     command = "A 1 AO 1 $v\n"
#     #     send_command(serial, command)
#     # end
# end

LibSerialPort.close(serial)
