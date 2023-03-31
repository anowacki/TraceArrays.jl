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

# Fallback methods
The following are defined for all `AbstractTraceArray`s, based on their containing
the above fields:
- `Seis.trace(::AbstractTraceArray)::AbstractMatrix`: The data matrix
- `Seis.nsamples(::AbstractTraceArray)::Int`: Number of samples
- `Seis.times(::AbstractTraceArray)`: Range giving times of each sample
- `Base.length(::AbstractTraceArray)::Int`: Number of channels
- `Base.getindex(::AbstractTraceArray, i::Int)::AbstractTrace`
- `Base.getindex(::AbstractTraceArray, indices)::AbstractVector{<:AbstractTrace}`

# Generated `Base` functions
`==` and `hash` are defined recursively on the fields of subtypes of
`AbstractTraceArray`.  These can be overloaded if needed.
"""
abstract type AbstractTraceArray <: AbstractTrace end

_geometry(::Event{T,P}) where {T,P} = P
_geometry(::Station{T,P}) where {T,P} = P

function Base.getindex(ta::AbstractTraceArray, i::Int)
    data = Seis.trace(ta)[:,i]
    t = Seis.Trace{typeof(ta.b),typeof(data),_geometry(ta.evt)}(
        ta.b, ta.delta, Seis.trace(ta)[:,i])
    t.evt = ta.evt
    t.sta = ta.sta[i]
    t.picks = ta.picks
    t.meta = ta.meta
    t
end
Base.getindex(ta::AbstractTraceArray, i) = [ta[ii] for ii in i]
Base.length(ta::AbstractTraceArray) = size(Seis.trace(ta), 2)

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

Seis.nsamples(ta::AbstractTraceArray) = size(Seis.trace(ta), 1)
Seis.trace(ta::AbstractTraceArray) = ta.data
