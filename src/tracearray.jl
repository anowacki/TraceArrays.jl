"""
    TraceArray{T,M,P} <: AbstractTraceArray

`TraceArray`s are sets of traces which are recorded on a strictly defined
array of seismic recorders.  This strict definition requires all channels
to have the same response and for all samples to be recorded at the same
time across array elements (meaning a uniform sampling interval as well).
For a single `TraceArray`, no channels can be partially present, meaning
all start and end on the same sample (although a `TraceArray` can be
constructed from channels where gaps have been filled).
"""
struct TraceArray{T,M,P} <: AbstractTraceArray
    b::T
    delta::T
    evt::Event{T,P}
    sta::Vector{Station{T,P}}
    data::M
    picks::Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}
    meta::Seis.SeisDict{Symbol,Any}
end

"""
    TraceArray(array_of_traces) -> ::TraceArray

Construct a `TraceArray` from an array of `AbstractTrace`s.

# Example
```
julia> using Seis

julia> TraceArray(filter(is_vertical, sample_data(:regional)))
```
"""
function TraceArray(t::AbstractArray{<:AbstractTrace})
    ntraces = length(t)
    if ntraces == 0
        throw(ArgumentError("cannot construct a TraceArray from an empty array"))
    end

    t1 = first(t)

    all(x -> Seis.nsamples(x) == Seis.nsamples(t1), t) ||
        throw(ArgumentError("all traces must be the same length"))
    all(x -> Seis.starttime(x) == Seis.starttime(t1), t) ||
        throw(ArgumentError("all traces must have the same start time"))
    all(x -> x.delta == t1.delta, t) ||
        throw(ArgumentError("all traces must have the same sampling interval"))
    # Use `===` to allow all to be missing
    all(x -> Seis.origin_time(x) === Seis.origin_time(t1), t) ||
        throw(ArgumentError("all traces must have the same origin time"))

    data = _matrix_from_traces(t)
end

"""
    _matrix_from_traces(array_of_traces) -> ::Matrix

Assuming that all the `AbstractTrace`s in `array_of_traces` have the same
number of samples, sampling rate and start time, construct a matrix of
"""
function _matrix_from_traces(t::AbstractArray{<:AbstractTrace})
    T = reduce(promote_type, typeof(Seis.starttime(tt)) for tt in t)
    M = Matrix{reduce(promote_type, eltype(tt) for tt in t)}
    P = reduce(promote_type, _geometry(tt.evt) for tt in t)

    t1 = first(t)

    # Has already been checked and not 0
    ntraces = length(t)
    # Already checked to be all the same
    npts = Seis.nsamples(t1)

    data = M(undef, npts, ntraces)
    for (matrix_col, tt) in enumerate(t)
        data[:,matrix_col] .= Seis.trace(tt)
    end

    picks = Seis.SeisDict(foldl(merge, t.picks))
    meta = Seis.SeisDict(foldl(merge, t.meta))

    TraceArray{T,M,P}(Seis.starttime(t1), t1.delta, t1.evt, t.sta, data, picks, meta)
end
