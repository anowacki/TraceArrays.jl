module TraceArraysMakieExt

import Seis
import TraceArrays

@static if isdefined(Base, :get_extension)
    import Makie
else
    import ..Makie
end
#=
    Record sections
=#
function Seis.plot_section!(
    ax::Makie.Axis,
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing,
    args...;
    kwargs...
)
    yvalues, _ = _get_y_values_and_label(t, y_values)
    Seis.plot_section!(ax, collect(t), yvalues, args...; kwargs...)
end
function Seis.plot_section!(
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing,
    args...;
    kwargs...
)
    yvalues, _ = _get_y_values_and_label(t, y_values)
    Seis.plot_section!(Makie.current_axis(), collect(t), yvalues, args...; kwargs...)
end

function Seis.plot_section(
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing,
    args...;
    axis=(),
    kwargs...
)
    yvalues, ylabel = _get_y_values_and_label(t, y_values)
    Seis.plot_section(collect(t), yvalues, args...; axis=(; ylabel, axis...), kwargs...)
end
function Seis.plot_section(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing,
    args...;
    axis=(),
    kwargs...
)
    yvalues, ylabel = _get_y_values_and_label(t, y_values)
    Seis.plot_section(gp, collect(t), yvalues, args...; axis=(; ylabel, axis...), kwargs...)
end

#=
    Heatmaps: time domain
=#
function TraceArrays.plot_heatmap!(
    ax::Makie.Axis,
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing;
    kwargs...
)
    yvalues, _ = _get_y_values_and_label(t, y_values)
    Makie.heatmap!(ax, Seis.times(t), yvalues, Seis.trace(t); kwargs...)
end

function TraceArrays.plot_heatmap(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing;
    axis=(),
    kwargs...
)
    yvalues, ylabel = _get_y_values_and_label(t, y_values)
    Makie.heatmap(gp, Seis.times(t), yvalues, Seis.trace(t);
        axis=(xlabel="Time / s", ylabel, axis...),
        kwargs...
    )
end

function TraceArrays.plot_heatmap(
    t::TraceArrays.AbstractTraceArray,
    y_values=nothing;
    figure=(),
    axis=(),
    kwargs...
)
    fig = Makie.Figure(; figure...)
    ax, hm = TraceArrays.plot_heatmap(fig[1,1], t, y_values; axis, kwargs...)
    Makie.FigureAxisPlot(fig, ax, hm)
end

#=
    Heatmaps: frequency-wavenumber domain
=#
function TraceArrays.plot_heatmap!(
    ax::Makie.Axis,
    t::TraceArrays.FourierDASArray;
    kwargs...
)
    Makie.heatmap!(ax, TraceArrays.wavenumbers(t), Seis.frequencies(t), abs.(Seis.trace(t))'; kwargs...)
end

function TraceArrays.plot_heatmap(
    gp::Union{Makie.GridPosition,Makie.GridSubposition},
    t::TraceArrays.FourierDASArray;
    axis=(),
    kwargs...
)
    Makie.heatmap(
        gp, TraceArrays.wavenumbers(t), Seis.frequencies(t), abs.(Seis.trace(t))';
        axis=(xlabel="Wavenumber / m⁻¹", ylabel="Frequency / Hz", axis...),
        kwargs...
    )
end

function TraceArrays.plot_heatmap(
    t::TraceArrays.FourierDASArray;
    figure=(),
    axis=(),
    kwargs...
)
    fig = Makie.Figure(; figure...)
    ax, hm = TraceArrays.plot_heatmap(fig[1,1], t; axis, kwargs...)
    Makie.FigureAxisPlot(fig, ax, hm)
end

# Add methods for Makie.heatmap[!] for time- and frequency-domain trace arrays
for AT in (TraceArrays.AbstractTraceArray, TraceArrays.AbstractFourierTraceArray)
    @eval begin
        Makie.heatmap(
            gp::Union{Makie.GridPosition,Makie.GridSubposition},
            t::$AT,
            args...;
            kwargs...
        ) = TraceArrays.plot_heatmap(gp, t, args...; kwargs...)
        Makie.heatmap(
            t::$AT,
            args...;
            kwargs...
        ) = TraceArrays.plot_heatmap(t, args...; kwargs...)
        Makie.heatmap!(
            ax::Makie.Axis,
            t::$AT,
            args...;
            kwargs...
        ) = TraceArrays.plot_heatmap!(ax, t, args...; kwargs...)
    end
end

#=
    Helpers
=#
"Return the y-coordinates at which to plot the traces and the axis label"
function _get_y_values_and_label(t::TraceArrays.AbstractTraceArray, y_values)
    if isnothing(y_values)
        _default_y_values_and_label(t)
    elseif y_values isa AbstractArray
        length(y_values) == length(t) ||
            throw(ArgumentError(
                "length of y values ($(length(y_values))) does not " *
                "match number of traces ($(length(t)))"
            ))
        y_values, ""
    elseif y_values == TraceArrays.distances
        y_values(t), "Distance / m"
    elseif y_values == Seis.distance_deg
        y_values.(t.evt, t.sta), "Epicentral distance / °"
    elseif y_values == Seis.distance_km
        y_values.(t.evt, t.sta), "Epicentral distance / km"
    elseif y_values isa Function
        y_values(t), ""
    elseif y_values isa Symbol
        getproperty.(t.meta, y_values), String(y_values)
    elseif y_values isa AbstractString
        if y_values == "index"
            1:length(t), "Trace index"
        else
            throw(ArgumentError("unrecognised y axis name '$y_values'"))
        end
    else
        throw(ArgumentError("unknown y_values value"))
    end
end

"The default y-coordinates and axis label for kinds of AbstractTraceArrays"
_default_y_values_and_label(t::TraceArrays.AbstractTraceArray) = 1:length(t), "Trace index"
function _default_y_values_and_label(t::TraceArrays.DASArray)
    ys = TraceArrays.distances(t)
    # Use km if range too large
    if abs(-(extrema(ys)...)) > 1000
        ys./1000, "Cable distance / km"
    else
        ys, "Cable distance / m"
    end
end

end # module
