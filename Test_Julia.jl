using GLMakie

positions = Observable(rand(Point2f, 10))
dragging = false
idx = 1

fig, ax, p = scatter(positions)

println(fig)
println(ax)
println(p)

on(events(fig).mousebutton, priority = 2) do event
    global dragging, idx
    if event.button == Mouse.left
        if event.action == Mouse.press
            plt, i = pick(fig)

            dragging = plt == p
            idx = i
            return Consume(dragging)
        elseif event.action == Mouse.release
            # Exit drag
            dragging = false
            return Consume(false)
        end
    end
    return Consume(false)
end

on(events(fig).mouseposition, priority = 2) do mp
    if dragging
        positions[][idx] = mouseposition(ax)
        notify(positions)
        return Consume(true)
    end
    return Consume(false)
end

fig