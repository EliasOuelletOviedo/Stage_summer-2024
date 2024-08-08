using GLMakie
using Observables

# Create a main window
main_fig = Figure(size = (600, 400))

# Function to create a pop-up window with a custom progress bar
function create_popup()
    # Create a new figure for the pop-up window
    popup_fig = Figure(size = (300, 150))

    # Create an axis to contain the progress bar
    ax = Axis(popup_fig[1, 1], title = "Progress", xticksvisible = false, yticksvisible = false)
    hidespines!(ax)
    hidedecorations!(ax)

    # Create a progress observable
    progress = Observable(0.0)

    # Initial polygon representing the progress bar (starting width = 0)
    progress_rect = @lift([Point2f(0, 0.4), Point2f($progress, 0.4), Point2f($progress, 0.6), Point2f(0, 0.6)])

    # Draw the polygon representing the progress bar
    poly!(ax, progress_rect, color = :green)

    # Add a label to show the progress percentage
    Label(popup_fig[2, 1], @lift("Progress: $(Int(round(100 * $progress)))%"))

    # Display the pop-up window
    display(popup_fig)

    # Function to run a task and update the progress bar
    function run_task()
        for i in 1:100
            sleep(0.1)  # Simulate work
            progress[] = i / 100
        end
    end

    # Run the task in a separate thread
    run_task()
end

# Add a button to the main window that triggers the pop-up
button = Button(main_fig[1, 1], label = "Run Task")
on(button.clicks) do _
    create_popup()
end

# Display the main window
display(main_fig)
