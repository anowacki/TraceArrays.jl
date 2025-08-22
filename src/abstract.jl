# Abstract types and their methods

"""
    AbstractTraceArray <: Seis.AbstractTrace

The abstract `AbstractTraceArray` type is the supertype for types which
correspond to sets of recordings of seismic data which are, strictly
interpreted, 'arrays'.  This means that all sensors have a common time
base and are of the same kind (i.e., have the same response).

The data underlying `AbstractTraceArray`s is usually stored as a matrix with
each column containing a continuous time series of data, and each row
representing the same time point at each array element (although the exact
implementation of `AbstractTraceArray` subtype need not be exactly like this).

# Indexing (`getindex`)
Indexing into a `AbstractTraceArray` with a single index gives you a single `Trace`
back; single-dimension slicing gives you a vector of `Trace`s.

!!! note
    In a future version of TraceArrays (or when it is merged into Seis.jl)
    slicing will instead return a new `AbstractTraceArray`.

# Required fields
Subtypes of `AbstractTraceArray` should have at least the following fields:
- `b`: Time in s relative to the origin time (if any) set in `evt.time`
- `delta`: Sampling interval in s
- `evt::Event`: Event relating to these data.  `.evt.time` is the time to which
  `.b` is relative.
- `sta::AbstractVector{<:Station}`: Set of stations, one for each channel
- `data::AbstractMatrix`: Data, where each column holds a single channel
- `picks::Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}`: Set of picks, where any
  one pick represents a single time across all channels.
- `meta::Seis.SeisDict{Symbol,Any}`: Metadata for all channels.

# Required methods
- `empty(::T) where {T<:AbstractTraceArray}`: Constructor for a new empty
  `AbstractTraceArray` where all fields are `deepcopy`ed apart from `.data`.
  This is intended to be used to implement copying functions (e.g., `Seis.decimate`)
  which do not need to first copy the whole `.data` field.

# Fallback methods
The following are defined for all `AbstractTraceArray`s, based on their containing
the above fields:
- `Seis.trace(::AbstractTraceArray)::AbstractMatrix`: The data matrix
- `Seis.nsamples(::AbstractTraceArray)::Int`: Number of samples
- `Seis.times(::AbstractTraceArray)`: Range giving times of each sample
- `Base.eltype(::AbstractTraceArray)::DataType`: Element type of data matrix
- `Base.length(::AbstractTraceArray)::Int`: Number of channels
- `Base.getindex(::AbstractTraceArray, i::Int)::AbstractTrace`
- `Base.getindex(::AbstractTraceArray, indices)::AbstractVector{<:AbstractTrace}`
- `Base.getindex(::AbstractTraceArray, ::Colon)::AbstractVector{<:AbstractTrace}`
- `Base.firstindex(::AbstractTraceArray)`
- `Base.lastindex(::AbstractTraceArray)`
- `(Base.empty(::T) where {T<:AbstractTraceArray})::T`
- `_check_abstracttracearray_consistency(::AbstractTraceArray, ::Vararg{AbstractTraceArray})`:
  Whether two or more arrays can be `append!`ed or `vcat`ted together.

# Generated `Base` functions
`==` and `hash` are defined recursively on the fields of subtypes of
`AbstractTraceArray`.  These can be overloaded if needed.
"""
abstract type AbstractTraceArray <: AbstractTrace end

_geometry(::Event{T,P}) where {T,P} = P
_geometry(::Station{T,P}) where {T,P} = P

function Base.getindex(ta::AbstractTraceArray, i::Integer)
    data = Seis.trace(ta)[:,begin+i-1]
    t = Seis.Trace{typeof(ta.b),typeof(data),_geometry(ta.evt)}(
        ta.b, ta.delta, data)
    t.evt = ta.evt
    t.sta = ta.sta[begin+i-1]
    t.picks = ta.picks
    t.meta = ta.meta
    t
end
Base.getindex(ta::AbstractTraceArray, i) = [ta[ii] for ii in i]
Base.getindex(ta::AbstractTraceArray, ::Colon) = [ta[i] for i in firstindex(ta):lastindex(ta)]
Base.axes(ta::AbstractTraceArray) = (Base.OneTo(length(ta)),)
Base.collect(ta::AbstractTraceArray) = [tt for tt in ta]
Base.CartesianIndices(ta::AbstractTraceArray) = CartesianIndices(axes(ta))
Base.LinearIndices(ta::AbstractTraceArray) = LinearIndices(axes(ta))
Base.length(ta::AbstractTraceArray) = size(Seis.trace(ta), 2)
Base.eltype(ta::AbstractTraceArray) = eltype(Seis.trace(ta))
Base.firstindex(ta::AbstractTraceArray) = first(LinearIndices(ta))
Base.lastindex(ta::AbstractTraceArray) = last(LinearIndices(ta))
Base.keys(ta::AbstractTraceArray) = LinearIndices(ta)
Base.iterate(ta::AbstractTraceArray, i=firstindex(ta)) = i > lastindex(ta) ? nothing : (ta[i], i + 1)
Base.pairs(ta::AbstractTraceArray) = Base.Pairs(values(ta), keys(ta))

@generated function Base.:(==)(t1::AbstractTraceArray, t2::AbstractTraceArray)
    fields = fieldnames(t1)
    types_equal = t1 == t2
    quote
        $types_equal || return false
        $([:(t1.$f == t2.$f || return false) for f in fields]...)
        true
    end
end

@generated function Base.hash(t::AbstractTraceArray, h::UInt)
    fields = fieldnames(t)
    quote
        $([:(h = hash(t.$f, h)) for f in fields]...)
        h
    end
end

"""
    append!(a::AbstractTraceArray, rest::AbstractTraceArray...; merge_meta=true) -> a

Append the traces in `rest` to `a`, requiring that they all have
consistent start times, sampling rates and trace lengths.

Metadata are merged together by default (unless `merge_meta` is `false`).
If any entries of the same name occur in both, then those in `a` are
used in preference; otherwise later entries in the argument list are kept.
"""
function Base.append!(a::AbstractTraceArray, rest::Vararg{AbstractTraceArray}; merge_meta=true)
    isempty(rest) && return a

    # Throw errors here if any problems
    _check_trace_consistency((a, rest...))

    a.data = hcat(Seis.trace(a), Seis.trace.(rest)...)
    for t in rest
        append!(a.sta, t.sta)
    end
    if merge_meta
        for t in rest
            a.meta = merge(t.meta, a.meta)
        end
    end
    a
end

"""
    empty(t::T) where {T<:AbstractTraceArray} -> t′::T

Create a copy of an `AbstractTraceArray` `t` where all fields apart
from the underlying trace data are copied to a new instance of `T`
using `Base.deepcopy`.

This function is not exported and is meant to be used by library code
to efficiently implement in- and out-of-place processing functions.

`t′.data` is set to be a 0×nchannels, matrix.  It is the caller's responsibility
to ensure that the matrix is set to be consistent with the other
internal fields (chiefly `.sta`) before `t′` is used further.
"""
@generated function Base.empty(t::T) where {T<:AbstractTraceArray}
    fields = fieldnames(t)
    quote
        $([f !== :data ? :($f = t.$f) : :($f = Matrix{eltype(t)}(undef, 0, 0)) for f in fields]...)
        # $(:(T($fields)))
        $(T)($([:($f) for f in fields]...))
    end
end

"""
    Base.vcat(a::AbstractTraceArray, rest::AbstractTraceArray...; merge_meta=true)

Append the traces in `rest` to `a`, requiring that they both have
consistent start times, sampling rates and trace lengths, returning
copies.

Metadata are merged together by default (unless `merge_meta` is `false`).
If any entries of the same name occur, then those in `a` are used in
preference; otherwise later entries in the argument list are kept.
"""
function Base.vcat(a::AbstractTraceArray, rest::Vararg{AbstractTraceArray}; merge_meta=true)
    isempty(rest) && return deepcopy(a)

    _check_trace_consistency((a, rest...))
    out = empty(a)
    out.data = hcat(Seis.trace(a), Seis.trace.(rest)...)
    out.sta = vcat(deepcopy(a.sta), map(x -> deepcopy(x.sta), rest)...)
    out.meta = deepcopy(a.meta)
    if merge_meta
        for t in rest
            out.meta = merge(deepcopy(t.meta), out.meta)
        end
    end
    out
end

"""
    _check_trace_consistency(ts::NTuple{N,AbstractTraceArray} where N)

Throw an error if the traces in `ts` cannot be sensibly joined together
into a new `AbstractTraceArray` with [`vcat`](@ref Base.vcat(::AbstractTraceArray, ::Vararg{AbstractTraceArray}))
or [`append!`](@ref Base.append(::AbstractTraceArray, ::Vararg{AbstractTraceArray})).

New subtypes of `AbstractTraceArray` may implement methods for this function
which impose additional checks on whether all the items in `ts` are
compatible with each other.
"""
function _check_trace_consistency(ts::NTuple{N,AbstractTraceArray} where N)
    _check_abstracttracearray_consistency(ts)
end

"""
    _check_abstracttracearray_consistency(ts)

Function which should throw if any `AbstractTraceArrays` cannot be
sensible merged together with `vcat` or `append`.
"""
function _check_abstracttracearray_consistency(ts)
    if length(ts) < 2
        throw(ArgumentError("number of trace arrays must be 2 or more"))
    end

    t1, rest = Iterators.peel(ts)

    all(x -> Seis.nsamples(x) == Seis.nsamples(t1), rest) ||
        throw(ArgumentError("all traces must be the same length"))
    all(x -> Seis.starttime(x) == Seis.starttime(t1), rest) ||
        throw(ArgumentError("all traces must have the same start time"))
    all(x -> x.delta == t1.delta, rest) ||
        throw(ArgumentError("all traces must have the same sampling interval"))
    # Use `===` to allow all to be missing
    all(x -> Seis.origin_time(x) === Seis.origin_time(t1), rest) ||
        throw(ArgumentError("all traces must have the same origin time"))
end

Seis.nsamples(ta::AbstractTraceArray) = size(Seis.trace(ta), 1)
Seis.trace(ta::AbstractTraceArray) = ta.data

"""
    AbstractFourierTraceArray <: AbstractFourierTrace

The abstract `AbstractFourierTraceArray` type is the supertype for
trace array types which are in the frequency (Fourier) domain.
"""
abstract type AbstractFourierTraceArray <: Seis.AbstractFourierTrace end
