"""
    read_optodas([T=Float64[, P=Seis.Geographic,]] file; kwargs...) -> ::DASArray

Read an HDF5 file written by an ASN OptoDAS interrogator.

Optionally specify the `DASArray` header floating point type `T`
(default `Float64`) and coordinate type (default `Seis.Geographic`).

# Keyword arguments
- `xdecimate = 1`: Decimate the number of channels when reading to every
  `xdecimate`th channel.  No spatial antialiasing is performed.
- `xlim`: A tuple of two values, giving the start and end distance (m) to read.
- `tlim`: A tuple of two values, giving the start and end time to read.  This
  may be either two real numbers in seconds, in which case the time is relative to
  the start of the trace, or two `Dates.AbstractDateTime`s such as
  `DateTime` giving the start and end time in UTC time.

!!! note
    Currently only files with a single region of interest are supported.
    If reading OptoDAS files with multiple regions of interest is required,
    please open an issue at the GtiHub repo for this package.
"""
function read_optodas(
    ::Type{T},
    ::Type{P},
    file::AbstractString;
    xdecimate::Integer=1,
    xlim::Union{Nothing, NTuple{2,Real}}=nothing,
    tlim::Union{Nothing, NTuple{2,Real}, NTuple{2,<:Dates.AbstractDateTime}}=nothing,
    header_only=false,
) where {T, P<:Seis.Position}
    # OptoDAS files don't use attributes; only groups and variables

    HDF5.h5open(file) do f
        file_version = f["fileVersion"][]
        if file_version != 8
            @warn "OptoDAS file version is $file_version; expected 8.  Errors may occur in reading."
        end

        header = f["header"]
        # Read all metadata, without touching the data itself
        meta = _optodas_meta(f)

        # Check units
        dimension_names = header["dimensionNames"][]
        if dimension_names != ["time", "distance"]
            @warn "OptoDAS data axes have dimensions of $(dimension_names), not [\"time\", \"distance\"]"
        end
        dimension_units = header["dimensionUnits"][]
        if dimension_units != ["s", "m"]
            @warn "OptoDAS data axes have units of $(dimension_units), not [\"s\", \"m\"]"
        end

        # Check data dimensions
        data_handle = f["data"]
        ndims(data_handle) == 2 ||
            throw(ArgumentError("`data` field of OptoDAS file has $(ndims(data_handle)) dimensions, not 2 as expected"))
        nchannels, npts = size(data_handle)
        if eltype(data_handle) != Float32
            @warn "OptoDAS data are not in Float32 precision.  Possible loss in conversion to Float32 upon reading."
        end

        demod_spec = f["demodSpec"]
        # ROI is 'region of interest': which section of the fibre to keep and
        # with what spacing.  Here we assume there is only one region of interest
        # and error if that is not the case.
        # TODO: Cope with multiple ROIs
        roi_start = demod_spec["roiStart"][]
        if length(roi_start) != 1
            throw(ArgumentError("file contains more than one ROI, which is not currently supported"))
        end

        # Can now assume a single ROI
        distances = f["cableSpec"]["sensorDistances"][]
        length(distances) == nchannels ||
            @warn "Unexpected number of distances $(length(distances)) compared to channels ($nchannels) in OptoDAS file"
        starting_distance = distances[begin]
        # Assumes evenly-spaced channels which increase in distance
        distance_spacing = (distances[end] - distances[begin])/(length(distances) - 1)

        time = _epochseconds_to_datetime(header["time"][])

        evt = Seis.Event{T,P}(; time)
        sta = [Seis.Station{T,P}() for _ in 1:nchannels]

        # Trace begin time will be zero relative to reference time
        # unless we are cutting; sorted out later
        b_zero = zero(T)
        delta = header["dt"][]

        picks = Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}()

        # Finally read in the data that we want, defining some
        # helper ranges first
        distances = starting_distance .+ (0:(nchannels - 1)).*distance_spacing
        times = b_zero .+ (0:(npts - 1)).*delta

        # Calculate trace start time; leave evt.time as is
        data, b = if header_only
            Matrix{Float32}(undef, 0, 0), b_zero
        elseif xdecimate == 1 && xlim == nothing && tlim == nothing
            permutedims(data_handle[]), b_zero
        else
            _xinds = _cut_indices(distances, xlim)
            xinds = _xinds[1]:xdecimate:_xinds[end]
            starting_distance = distances[xinds[1]]

            # Only approximate cutting with dates if sampling rate not a whole
            # number of milliseconds
            tinds = if tlim isa NTuple{2, Dates.AbstractDateTime}
                dates = evt.time .+ Dates.Millisecond.(floor.(Int, times.*1000))
                _cut_indices(dates, tlim)
            else
                _cut_indices(times, tlim)
            end

            b = times[tinds[1]]

            permutedims(data_handle[xinds,tinds]), b
        end

        TraceArrays.DASArray{T,Matrix{Float32},P}(
            b, delta, starting_distance, distance_spacing, evt, sta,
            data, picks, meta
        )
    end
end
read_optodas(::Type{T}, file::AbstractString; kwargs...) where T = read_optodas(T, Seis.Geographic{T}, file; kwargs...)
read_optodas(file::AbstractString; kwargs...) = read_optodas(Float64, file; kwargs...)

"""
Find the indices of `arr` whose values are between `lim[1]` and `lim[2]`.
This strictly only inludes elements of `arr` greater than or equal to `lim[1]`
and less than or equal to `lim[2]`.
"""
_cut_indices(arr, lim::Nothing) = eachindex(arr)
function _cut_indices(arr, lim)
    i1 = max(searchsortedfirst(arr, lim[1]), 1)
    i2 = min(searchsortedlast(arr, lim[2]), length(arr))
    i1:i2
end

"""
Convert fractional seconds since Unix epoch to proper `Dates.DateTime`.
Assumes this is always a whole number of milliseconds and warn if not.
"""
function _epochseconds_to_datetime(seconds::Real)
    if rem(seconds*1000, 1) != 0
        @warn "Unexpected non-integer number of milliseconds since Unix epoch in OptoDAS header"
    end
    ms = round(Int, (seconds%1)*1000)
    Dates.unix2datetime(floor(Int, seconds)) + Dates.Millisecond(ms)
end

function _optodas_meta(f)
    # Fill common DASArray meta information
    header = f["header"]
    meta = Seis.SeisDict{Symbol,Any}()
    meta.gauge_length = header["gaugeLength"][]

    # Save all OptoDAS information
    meta.optodas = Seis.SeisDict{Symbol,Any}()
    _fill_meta!(meta.optodas, f)

    meta
end

function _fill_meta!(dict::Seis.SeisDict, h5)
    for key in keys(h5)
        key == "data" && continue
        symbol_key = Symbol(key)

        val = h5[key]

        if typeof(val) == HDF5.Group
            dict[symbol_key] = Seis.SeisDict{Symbol,Any}()
            _fill_meta!(dict[symbol_key], val)
        else
            dict[symbol_key] = val[]
        end
    end
end
