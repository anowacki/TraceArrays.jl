function Seis.decimate!(t::AbstractTraceArray, n::Integer; antialias=true)
    n == 1 && return t

    newdata = _decimate_core(t, n, antialias)
    t.data = newdata
    t.delta *= n
    t
end

function Seis.decimate(t::AbstractTraceArray, n::Integer; antialias=true)
    n == 1 && return deepcopy(t)

    newdata = _decimate_core(t, n, antialias)
    t′ = empty(t)
    t′.data = newdata
    t′.delta *= n
    t′
end

function _decimate_core(t::AbstractTraceArray, n, antialias)
    1 <= n || throw(ArgumentError("n must be greater than 0 (supplied $n)"))
    olddata = Seis.trace(t)

    newdata = if antialias
        DSP.resample(olddata, 1//n; dims=1)
    else
        olddata[begin:n:end,:]
    end
    newdata
end

# FIXME: Implement filtering in a more general and efficient way
function Seis.bandpass!(t::AbstractTraceArray, f1, f2;
    poles=2, twopass=false, kind=DSP.Butterworth(poles)
)
    data = Seis.trace(t)
    T = eltype(data)
    temp = Seis.Trace(t.b, t.delta, Vector{T}(undef, Seis.nsamples(t)))
    for icol in axes(Seis.trace(t), 2)
        Seis.trace(temp) .= view(data, :, icol)
        Seis.bandpass!(temp, f1, f2; poles, twopass, kind)
        data[:,icol] .= Seis.trace(temp)
    end
    t
end

function Seis.lowpass!(t::AbstractTraceArray, f;
    poles=2, twopass=false, kind=DSP.Butterworth(poles)
)
    data = Seis.trace(t)
    T = eltype(data)
    temp = Seis.Trace(t.b, t.delta, Vector{T}(undef, Seis.nsamples(t)))
    for icol in axes(Seis.trace(t), 2)
        Seis.trace(temp) .= view(data, :, icol)
        Seis.lowpass!(temp, f; poles, twopass, kind)
        data[:,icol] .= Seis.trace(temp)
    end
    t
end

function Seis.highpass!(t::AbstractTraceArray, f;
    poles=2, twopass=false, kind=DSP.Butterworth(poles)
)
    data = Seis.trace(t)
    T = eltype(data)
    temp = Seis.Trace(t.b, t.delta, Vector{T}(undef, Seis.nsamples(t)))
    for icol in axes(Seis.trace(t), 2)
        Seis.trace(temp) .= view(data, :, icol)
        Seis.highpass!(temp, f; poles, twopass, kind)
        data[:,icol] .= Seis.trace(temp)
    end
    t
end

function Seis.resample!(t::AbstractTraceArray; n=nothing, delta=nothing)
    rate = _resample_rate(t, n, delta)
    rate == 1 && return t # Nothing to do
    t.data = DSP.resample(Seis.trace(t), rate; dims=1)
    t.delta /= rate
    t
end

function Seis.resample(t::AbstractTraceArray; n=nothing, delta=nothing)
    rate = _resample_rate(t, n, delta)
    rate == 1 && return deepcopy(t) # Nothing to do
    t′ = empty(t)
    t′.data = DSP.resample(Seis.trace(t), rate; dims=1)
    t′.delta /= rate
    t′
end

function _resample_rate(t::AbstractTraceArray, n, delta)
    if (delta === nothing && n === nothing) || (delta !== nothing && n !== nothing)
        throw(ArgumentError("one and only of `delta` and `n` must be given"))
    end
    rate = if delta !== nothing
        t.delta/delta
    else # n !== nothing
        n
    end
end

"""
    remove_median_trace!(t::AbstractTraceArray) -> t

Remove the median trace from all channels in `t`.

The operation finds the median value at each time sample across all channels
and subtracts this from all channels, modifying the original trace.

See also [`remove_median_trace`](@ref).
"""
function remove_median_trace!(t::AbstractTraceArray)
    data = Seis.trace(t)
    for ichannel in axes(data, 1)
        data[ichannel,:] .-= median(@view(data[ichannel,:]))
    end
    t
end

"""
    remove_median_trace(t::AbstractTraceArray) -> t′

Out-of-place version of [`remove_median_trace!`](@ref).
"""
remove_median_trace(t::AbstractTraceArray) = remove_median_trace!(deepcopy(t))
