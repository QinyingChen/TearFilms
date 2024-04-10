gx, gy, sol, h, p, c = twodim_solve(50, 50, (0, 3));

##
time = Observable(0.0)

cr = extrema([h(3);h(0)])
fig, ax, plt = heatmap(gx.t, gy.t, @lift(h($time)), colorrange=cr, interpolate=true)
ax.autolimitaspect[] = 1
Colorbar(fig[1,2], plt)

framerate = 30
timestamps = range(0, 3, step=1/framerate)

record(fig, "tmp.mp4", timestamps;
framerate = framerate) do t
    time[] = t
    ax.title[] = "t = $(round(t, digits = 1))"
end
