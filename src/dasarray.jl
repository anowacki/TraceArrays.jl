"""
    DASArray{T,M,P<:Seis.Position} <: AbstractTraceArray

`AbstractTraceArray` holding DAS data.

This type differs from `TraceArray` in that it records the information on
gauge (channel) spacing and distance, which is useful for some things.

---

    DASArray(; kwargs...)

Return a new, empty object for DAS data.  Individual fields may be defined using
keyword arguments.

!!! note
    The keyword argument constructor for `DASArray`s is not a stable part of the
    API and may change.

# Keyword arguments
- `data`: `AbstractMatrix` of data, where columns contain continuous
  evenly-sampled recordings at each distance, and each column is recorded
  a constant distance from the last.  These columns therefore correspond to
  the gauges used in DAS recordings.  The number of channels (columns) is used
  to create the correct number of stations in the `.sta` field.
- `starting_distace`: Distance from the interrogator of the first channel (m)
- `distance_spacing`: Spacing between channels (m)
- `b`: Start time (relative to `.evt.time`) of the first sample (s)
- `delta`: Sampling interval (s)
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
    data=Array{T}(undef, 0, 0),
    b=0,
    delta=1,
    starting_distance=0,
    distance_spacing=1,
    evt=Event{T,P}(),
    sta=[Station{T,P}() for _ in axes(data, 2)],
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
    DASArray(::AbstractTraceArray, starting_distance, distance_spacing) -> ::DASArray

Construct a `DASArray` from another kind of `AbstractTraceArray` (like a
`TraceArray`).  Supply the distance of the first channel in m as
`starting_distance` and the channel spacing in `distance_spacing`, also
in m.
"""
function DASArray(t::AbstractTraceArray, starting_distance, distance_spacing)
    M = typeof(Seis.trace(t))
    T = typeof(Seis.starttime(t))
    P = _geometry(t.evt)
    DASArray{T,M,P}(Seis.starttime(t), t.delta, starting_distance, distance_spacing,
        t.evt, t.sta, Seis.trace(t), t.picks, t.meta)
end

# Specialised constructor for `DASArray`s since `DASArray <: AbstractTraceArray`
DASArray(t::DASArray, starting_distance, distance_spacing) = t

"""
    DASArray(array_of_traces::AbstractArray{<:AbstractTrace}, starting_distance, distance_spacing)

Construct a `DASArray` from an array of single-channel traces, providing the
distance to the first channel in m as `starting_distance`, and the channel
spacing in m as `distance_spacing`.
"""
function DASArray(ts::AbstractArray{<:AbstractTrace}, starting_distance, distance_spacing)
    DASArray(TraceArray(ts), starting_distance, distance_spacing)
end

"""
    distances(da::DASArray) -> dists

Return a range giving the distances in m along the fibre of each channel in `da`.
"""
distances(t::DASArray) = t.starting_distance .+ (0:(size(Seis.trace(t), 2) - 1)).*t.distance_spacing

"""
    read_febus(file; kwargs...) -> ::DASArray{Float64, Matrix{Float32}, Seis.Geographic{Float64}}

Read DAS data from `file` in the Febus HDF5 format.

# Keyword arguments
- `xlim=(-Inf, Inf)`: Read only a range of distances (m)
- `xdecimate=1`: Decimate the number of channels read in.
- `tlim=($(typemin(DateTime)), $(typemax(DateTime)))`: Read data between two `DateTime`s
- `blocks=(1,nblocks)`: Read in only data blocks within the given range.
  (Data blocks are typically 1 s long, but can be any length.)
"""
function read_febus(file; kwargs...)
    data = FebusTools.read_hdf5(file; kwargs...)
    DASArray(data)
end

"""
    cut_distance!(t::DASArray, x1, x2; warn=true, allowempty=false) -> t

Cut out channels in the `DASArray` `t` which lie between distances `x1` m
and `x2` m.  Only channels which are at exactly `x1` m or more, and exactly
`x2` m or less are retained.

If `allowempty` is `true`, then an empty `DASArray` can be returned when
`x1` and `x2` are either both below or above the cable distance range.
Otherwise, an error is thrown.

By default, warnings are printed if `x1` or `x2` lie outside the existing
range of distances for `t` and the start or end distance is used instead
to perform the cut.  Use `warn=false` to turn these warnings off.

See also: [`cut_distance`](@ref).
"""
function cut_distance!(t::DASArray, x1, x2; warn=true, allowempty=false)
    ib, ie = _cut_distance_indices(t, x1, x2; warn, allowempty)
    t.data = t.data[:,ib:ie]
    old_distances = distances(t)
    t.starting_distance = ib in eachindex(old_distances) ? old_distances[ib] : x1
    t.sta = t.sta[ib:ie]
    t
end

"""
    cut_distance(t::DASArray, x1, x2; warn=true, allowempty=false) -> t′

Out-of-place version of [`cut_distance!`](@ref).
"""
function cut_distance(t::DASArray{T,M,P}, x1, x2; warn=true, allowempty=false) where {T,M,P}
    ib, ie = _cut_distance_indices(t, x1, x2; warn, allowempty)
    old_distances = distances(t)
    starting_distance = ib in eachindex(old_distances) ? old_distances[ib] : x1
    sta = t.sta[ib:ie]
    data = Seis.trace(t)[:,ib:ie]
    DASArray(t; starting_distance, sta, data)
end

function _cut_distance_indices(t::DASArray, b, e; warn=true, allowempty=false)
    start_dist = first(distances(t))
    end_dist = last(distances(t))
    (b === missing || e === missing) && throw(ArgumentError("Start or end cut distance is `missing`"))
    e < b && throw(ArgumentError("End cut distance ($e m) is before starting cut ($b m)"))
    if b > end_dist || e < start_dist
        if !allowempty
            b > end_dist &&
                throw(ArgumentError("Beginning cut distance $b m is beyond end of data ($(end_dist) m)."))
            e < start_dist &&
                throw(ArgumentError("End cut distance $e m is earlier than start of data ($(end_dist) m)."))
        end
        empty!(t.t)
        t.b = b
        return t
    end
    if b < start_dist
        warn && @warn("Beginning cut distance $b m is before start of data.  Setting to $(start_dist) m.")
        b = start_dist
    end
    if e > end_dist
        warn && @warn("End cut distance $e m is after end of data.  Setting to $(end_dist) m.")
        e = end_dist
    end
    ib = round(Int, (b - start_dist)/t.distance_spacing) + 1
    ie = length(t) - round(Int, (end_dist - e)/t.distance_spacing)
    ib, ie
end


function decimate_distance!(t::DASArray, n::Integer)
    t.data = t.data[:,begin:n:end]
    t.distance_spacing *= n
    t.sta = t.sta[begin:n:end]
    t
end

function decimate_distance(t::DASArray, n::Integer)
    sta = t.sta[begin:n:end]
    distance_spacing = t.distance_spacing*n
    data = Seis.trace(t)[:,begin:n:end]
    DASArray(t; distance_spacing, sta, data)
end

"""
    integrate_distance(t::DASArray, n=3) -> t′

Integrate the data in `t` with respect to distance along the
cable.  If data are in units of strain, therefore, they are
converted to units of displacement (and similarly for strain rate
to velocity).  No unit conversion is made, so data which are in
nanostrain convert to nm, and likewise nanostrain/s becomes
nm/s.

`n` determines how many points are used in the integration and
must be an odd number.  This determines the number of points
used for the trapezoidal integration and must be an odd number.

The number of channels is reduced by `n - 1`.  For example,
when `n` is 3, `t′` has 2 fewer channels than before, and the
first and last channels in `t` are removed from `t′`.

# Example
```
julia> t = DASArray(; data=[0. 1 0 -1; 1 2 -2 1], b=0, delta=1, starting_distance=0, distance_spacing=0.5);

julia> integrate_distance(t)
2×2 Matrix{Float64}:
 0.5    0.0
 0.75  -0.25
```
"""
function integrate_distance(t::DASArray, n::Integer=3)
    isodd(n) && n > 1 || throw(ArgumentError("`n` must be and odd number 3 or more"))

    nchannels = length(t)
    nchannels′ = nchannels - n + 1
    npts = Seis.nsamples(t)
    t′ = empty(t)
    data = Seis.trace(t)
    data′ = similar(data, npts, nchannels′)
    t′.data = data′

    u = similar(data, nchannels)
    U = similar(data, nchannels′)

    # Saves a divide in each timestep
    spacing_by_2 = t.distance_spacing/2

    for i in axes(data, 1)
        u = data[i,:]
        U = data′[i,:]
        _integrate!(U, u, spacing_by_2, n)
        data′[i,:] .= U
    end

    # Update station information
    t′.sta = t′.sta[(begin + n÷2):(end - n÷2)]
    t′.starting_distance += (n÷2)*t′.distance_spacing

    t′
end

"""
    _integrate!(U, u, spacing_by_2, n) -> U

Integrate the values of `u` using an `n`-point trapezium rule, putting
the output into `U`.  `U` must be the correct size, i.e., it must have
`n - 1` points fewer than `u`.  `spacing_by_2` is half the even spacing
between points in `u`.

# Example
```
julia> u = 1:7;

julia> U = similar(u, 5);

julia> TraceArrays._integrate!(U, u, 10/2, 3)
5-element Vector{Float64}:
  40.0
  60.0
  80.0
 100.0
 120.0
```
"""
function _integrate!(U, u, spacing_by_2, n)
    # Special cases which are quicker
    if n == 3
        for i in eachindex(U)
            U[i] = spacing_by_2*(u[i] + 2*u[i+1] + u[i+2])
        end

    elseif n == 5
        for i in eachindex(U)
            U[i] = spacing_by_2*(u[i] + 2*(u[i+1] + u[i+2] + u[i+3]) + u[i+4])
        end

    # General case: much slower, though `@inbounds` speeds up by 35%
    else
        for i in eachindex(U)
            # First point
            U[i] = u[i]
            # Inner points.  Only runs twice for n == 4.
            for j in (i + 1):(i + n - 2)
                U[i] += 2*u[j]
            end
            # Last point
            U[i] += u[i+n-1]
            U[i] *= spacing_by_2
        end
    end
    U
end
