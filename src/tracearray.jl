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
mutable struct TraceArray{T,M,P} <: AbstractTraceArray
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
    if isempty(t)
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

    data, T, M, P = _matrix_from_traces(t)
    picks = Seis.SeisDict(foldl(merge, t.picks))
    meta = Seis.SeisDict(foldl(merge, t.meta))

    TraceArray{T,M,P}(Seis.starttime(t1), t1.delta, t1.evt, t.sta, data, picks, meta)
end

"""
    _matrix_from_traces(array_of_traces) -> ::Matrix, T, M, P

Assuming that all the `AbstractTrace`s in `array_of_traces` have the same
number of samples, sampling rate and start time, construct a matrix of
the data where each column contains the samples from the corresponding
trace.

Also return `T`, `M` and `P`, the type parameters for the `TraceArray`
which can be formed with the matrix.
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

    data, T, M, P
end
