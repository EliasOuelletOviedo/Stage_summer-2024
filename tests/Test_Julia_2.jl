using GLMakie
using Observables

# Create observables for data points
x_data = Observable(Float64[])
y_data = Observable(Float64[])

# Create an observable for the running status
running = Observable(false)

# Define a function that generates new points
function generate_points()
    t = 0.0
    while running[]
        push!(x_data[], t)
        push!(y_data[], sin(t))
        notify(x_data)
        notify(y_data)
        t += 0.1
        sleep(0.1)
    end
end

# Create the figure and plot
fig = Figure(size = (1000, 750))
ax = Axis(fig[1, 1:4], title="Real-time plot")
line = lines!(ax, x_data, y_data)

# Create start and stop buttons
button_start = Button(fig[2, 2], label="Start")
button_stop = Button(fig[2, 3], label="Stop")

# Start the function when the start button is pressed
on(button_start.clicks) do _
    running[] = false
    sleep(0.1)
    running[] = true
    @async generate_points()  # Start a new task
end

# Stop the function when the stop button is pressed
on(button_stop.clicks) do _
    running[] = false
end

fig
