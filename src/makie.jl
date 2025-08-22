"""
    plot_heatmap([gridposition,] t::DASArray; figure=(), axis(), kwargs...)
    plot_heatmap([gridposition,] t::FourierDASArray; figure=(), axis(), kwargs...)

Plot the trace values of a `DASArray`, or the absolute values of
a `FourierDASArray`'s coefficients, in terms of frequency and wavenumber
on via `Makie.heatmap`.

If the first argument is a grid position (`Makie.GridPosition` or
`Makie.GridSubposition`), then a new axis is created at that point within
the layout.  Keyword arguments to the `Makie.Axis` constructor can be
passed as a named tuple in `axis`.  This method returns a `Makie.AxisPlot`
which can be iterated to retrieve the axis handle and plot handle.

Without a grid position, a new figure
is created, in which case the `figure` keyword argument can be used to
pass a named tuple of keyword arguments to the `Makie.Figure` constructor.
This method returns a `Makie.FigureAxisPlot` containing handles to all
three of the figure, axis and plot object.

Remaining keyword arguments `kwargs` are passed to the call to
`Makie.heatmap`.

This function can only be used after loading a
[Makie](https://docs.makie.org/stable/)
[backend](https://docs.makie.org/stable/explanations/backends/backends),
e.g. by doing `using GLMakie`, `using CairoMakie`, and so on.
"""
function plot_heatmap end

"""
    plot_heatmap!(axis::Makie.Axis, t::FourierDASArray; kwargs...)

Plot the trace values of a `DASArray`, or the absolute values of
a `FourierDASArray`'s coefficients, in terms of frequency and wavenumber
on via `Makie.heatmap`.  The plot is added to an existing `Makie.Axis`
`axis`.  `kwargs` are passed to `Makie.heatmap`.

This function can only be used after loading a
[Makie](https://docs.makie.org/stable/)
[backend](https://docs.makie.org/stable/explanations/backends/backends),
e.g. by doing `using GLMakie`, `using CairoMakie`, and so on.
"""
function plot_heatmap! end
