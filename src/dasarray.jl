"""
    DASArray{T,M,P<:Seis.Position} <: AbstractTraceArray

`AbstractTraceArray` holding DAS data.

This type differs from `TraceArray` in that it records the information on
gauge (channel) spacing and distance, which is useful for some things.

---

    DASArray()

Return a new, empty object for DAS data.  Individual fields may be defined using
keyword arguments.

!!! note
    The keyword argument constructor for `DASArray`s is not a stable part of the
    API and may change.
"""
mutable struct DASArray{T,M,P} <: AbstractTraceArray
    b::T
    delta::T
    starting_distance::T
    distance_spacing::T
    evt::Event{T,P}
    sta::Vector{Station{T,P}}
    data::M
    picks::Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}
    meta::Seis.SeisDict{Symbol,Any}
end

function DASArray(;
    T=Float64,
    P=Seis.Geographic{T},
    b=0,
    delta=1,
    starting_distance=0,
    distance_spacing=1,
    evt=Event{T,P}(),
    sta=Station{T,P}[],
    data=Array{T}(undef, 0, 0),
    picks=Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}(),
    meta=Seis.SeisDict{Symbol,Any}(),
    M=typeof(data),
)
    DASArray{T,M,P}(b, delta, starting_distance, distance_spacing, evt, sta, data, picks, meta)
end

"""
    DASArray(t::DASArray{T,M,P}; kwargs...) -> t′

Create a new `DASArray` from an existing one, but where the trace data
is replaced by `data`.  The data matrix type (`M` parameter in the
`DASArray` type) is replaced by the type of `data`.

If the number of columns in `data` is not the same as the number of channels
in `t`, then the vector of stations in `t′.sta` is adjusted to be the same
length as `data`.  If there are more channels than before, the `t′.sta` vector
will contain `#undef` elements where the additional channels are present.
Accessing these elements before properly setting them will throw an
`UndefRefError`.
"""
function DASArray(t::DASArray{T0,M0,P0};
    T=T0,
    P=P0,
    b=Seis.starttime(t),
    delta=t.delta,
    starting_distance=first(distances(t)),
    distance_spacing=step(distances(t)),
    evt=deepcopy(t.evt),
    data=deepcopy(Seis.trace(t)),
    sta=resize!(deepcopy(t.sta), size(data, 2)),
    picks=deepcopy(t.picks),
    meta=deepcopy(t.meta),
    M=typeof(data)
) where {T0,M0,P0}
    DASArray{T,M,P}(b, delta, starting_distance, distance_spacing, evt, sta, data, picks, meta)
end

"""
    DASArray(t::FebusTools.FebusData)

Construct a `DASArray` from data read from a Febus A1 interrogator.

Instrument metadata is placed in the trace's `.meta.febus` field.
Note that some metadata items (such as `:Extent`) are not modified from
the file read, and so may not be in sync with the data returned.
"""
function DASArray(t::FebusTools.FebusData)
    data = t.data
    b = t.times[begin]
    delta = t.metadata[:Spacing][2]/1000 # Time sample spacing in ms
    evt = Seis.Event(; time=first(t.dates))
    starting_distance = first(t.distances)
    distance_spacing = step(t.distances)
    stas = [Seis.Station() for _ in axes(data, 2)]
    picks = Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{Float64}}()
    meta = Seis.SeisDict{Symbol,Any}()
    meta.febus = Seis.SeisDict{Symbol,Any}(t.metadata)
    DASArray(b, delta, starting_distance, distance_spacing, evt, stas, data, picks, meta)
end

"""
    distances(da::DASArray) -> dists

Return a range giving the distances in m along the fibre of each channel in `da`.
"""
distances(t::DASArray) = t.starting_distance .+ (0:(size(Seis.trace(t), 2) - 1)).*t.distance_spacing

"""
    read_febus(file; kwargs...) -> ::DASArray{Float64, Matrix{Float64}, Seis.Geographic{Float64}}

Read DAS data from `file` in the Febus HDF5 format.

# Keyword arguments
- `xlim=(-Inf, Inf)`: Read only a range of distances (m)
- `tlim=($(typemin(DateTime)), $(typemax(DateTime)))`: Read data between two `DateTime`s
- `blocks=(1,nblocks)`: Read in only data blocks within the given range.
  (Data blocks are typically 1 s long, but can be any length.)
"""
function read_febus(file; kwargs...)
    data = FebusTools.read_hdf5(file; kwargs...)
    DASArray(data)
end
