"""
    FourierDASArray <: AbstractFourierTraceArray

The Fourier transform of a [`DASArray`](@ref).

Because `DASArray`s by definition are equally sampled in space along
the fibre, `FourierDASArrays` represent the 2D Fourier transform of
the data sampled.  `FourierDASArray`s contain all the information needed
on the original number of samples (in time) and channels (in space).
The fields `b`, `delta` (frequency spacing in Hz), `starting_distance`
`distance_spacing`, `sta`, `evt`, `picks` and `meta` remain part of
the public API of the type.

The usual way to create `FourierDASArray`s is by calling [`fft`](@ref)
on a `DASArray`.  To convert back to a `DASArray`, call [`ifft`](@ref)

!!! note
    `FourierDASArray`s share many fields with the `DASArray` they came from, including
    the station and event, plus picks and metadata.  These fields are **not**
    copied across, and are instead references to the same `Station`, `Event`,
    and so on.  Therefore, any changes to the `FourierDASArray` will be reflected
    in the corresponding `DASArray`.  To avoid this, do `f = fft(deepcopy(t))`.

---
See also: [`DASArray`](@ref), [`fft`](@ref), [`ifft`](@ref).
"""
mutable struct FourierDASArray{
    T<:AbstractFloat,
    M<:AbstractMatrix{<:Complex},
    P<:Seis.Position{T}
} <: AbstractFourierTraceArray
    b::T
    # Frequency spacing in Hz
    delta::T
    starting_distance::T
    distance_spacing::T
    nsamples::Int64
    sta::Vector{Station{T,P}}
    evt::Event{T,P}
    data::M
    picks::Seis.SeisDict{Union{Int,Symbol}, Seis.Pick{T}}
    meta::Seis.SeisDict{Symbol,Any}

    function FourierDASArray{T,M,P}(
        b, delta, starting_distance, distance_spacing, nsamples,
        sta, evt, data, picks, meta
    ) where {T,M,P}
        delta > 0 || throw(ArgumentError("delta cannot <= 0"))
        nsamples >= 0 || throw(ArgumentError("nsamples cannot be negative"))
        new{T,M,P}(
            b, delta, starting_distance, distance_spacing, nsamples,
            sta, evt, data, picks, meta
        )
    end
end

"""
    FourierDASArray(; b, delta, data, nsamples=(2*length(data) - 1), starting_distance, distance_spacing, evt=Event(), sta=Station(), picks=nothing, meta=nothing)

Create a `FourierDASArray` directly from a set of Fourier coefficients.
Users will usually construct a `FourierDASArray` by calling [`fft`](@ref)
on a `Trace` instead of using this constructor.

# Keyword arguments
- `b`: The starting time in s of the original recording this frequency
  domain trace represents.
- `delta`: The **original sampling interval in s** of the equivalent time
  domain trace.
- `data`: An `AbstractMatrix{<:Complex}` containing the set of Fourier
  coefficients for this frequency domain trace.  Note that this is a
  'one-sided' set, where the first index of the first dimension corresponds
  to 0 Hz and the final index is the Nyquist frequency.  This is because
  `DASArray`s represent real quantities.  In the second dimension, the first
  index corresponds to the most negative wavenumber, the central index a
  wavenumber of 0, and the final index the largest wavenumber, where the
  most negative and largest wavenumbers have the same absolute value.
- `nsamples`: The number of time domain samples in the origin time
  domain trace which this frequency domain trace represents.  This allows
  the `FourierDASArray` to be converted back to the original `DASArray` with no
  loss of the number of points.
- `starting_distance`: Distance along cable of the first channel (m).
- `distance_spacing`: Spacing between channels along the cable (m).
- `evt`: An `Event` which defines the source and origin time for the data.
  Normally this should be taken from the `Trace` being used to construct
  the `FourierDASArray`.
- `sta`: A `Vector{<:Station}` which defines the recording stations for the data.
  Normally this should be taken from the `DASArray` being used to construct
  the `FourierDASArray`.
"""
function FourierDASArray{T,M,P}(;
    b, delta, data,
    nsamples=(2*size(data, 1) - 1),
    starting_distance,
    distance_spacing,
    evt=Event{T,P}(),
    sta=[Station{T,P}() for _ in 1:(size(data, 2)÷2 + 1)],
    picks=Seis.SeisDict{Union{Symbol,Int},Seis.Pick{T}}(),
    meta=Seis.SeisDict{Symbol,Any}()
) where {T,M,P}
    FourierDASArray{T,M,P}(
        b, 1/(nsamples*delta), starting_distance, distance_spacing, nsamples,
        sta, evt, data, picks, meta
    )
end

FourierDASArray{T,M}(; kwargs...) where {T,M} = FourierDASArray{T,M,Seis.Geographic{T}}(; kwargs...)
FourierDASArray{T}(; kwargs...) where {T} =
    FourierDASArray{T,Matrix{Complex{T}},Seis.Geographic{T}}(; kwargs...)
FourierDASArray(; kwargs...) =
    FourierDASArray{Float64,Matrix{Complex{Float64}},Seis.Geographic{Float64}}(; kwargs...)

@eval begin
    function Base.hash(f::FourierDASArray, h::UInt)
        $([:(h = hash(f.$f, h)) for f in fieldnames(FourierDASArray)]...)
        h
    end

    function Base.:(==)(a::FourierDASArray, b::FourierDASArray)
        all(isequal(getfield(a, f), getfield(b, f)) for f in $(fieldnames(FourierDASArray)))
    end
end

Base.broadcastable(f::FourierDASArray) = Ref(f)

distances(f::FourierDASArray) = f.starting_distance .+ (0:(size(Seis.trace(f), 2) - 1)).*f.distance_spacing
Seis.nfrequencies(f::FourierDASArray) = size(trace(f), 1)

"""
    Seis.fft(t::DASArray)

Convert the DAS data `t` into its equivalent frequency domain data `f` by
performing a Fourier transform.  The object returned is a [`FourierDASArray`](@ref).

# Example
```
julia> v_app = 5 # Apparent velocity m/s;

julia> freq = 0.5 # Wave frequency, Hz;

julia> delta = 0.05 # Sampling interval;

julia> distance_spacing = 0.05 # Channel spacing, m;

julia> t = DASArray(;
           b=0, delta, starting_distance=0, distance_spacing,
           data=sin.(2π*freq.*((0:201).*delta .- (0:201)'.*distance_spacing./v_app))
       ); # Sine wave of frequency 1 Hz with apparent velocity of 5 m/s;

julia> T = fft(t);

julia> indmax = argmax(abs.(trace(T))) # Get index of maximum power
CartesianIndex(6, 101)

julia> f = frequencies(T)[indmax[1]] # Maximum frequency in Hz is ~1 as expected
0.495049504950495

julia> wavenumbers(T) # Wavenumbers corresponding to second dimension of `trace(f)`
-9.999999999999998:0.099009900990099:9.900990099009901

julia> λ⁻¹ = -wavenumbers(T)[indmax[2]] # Wavenumber of wave in m⁻¹ is -1/-5 ≈ 0.2 as expected
0.09900990099009732

julia> f/λ⁻¹ # Wave velocity is frequency over wavenumber; ~5 m/s as expected
5.000000000000084
```

See also: [`ifft`](@ref).
"""
function Seis.fft(t::DASArray{T,<:AbstractMatrix{TV},P}) where {T,TV,P}
    data = FFTW.fftshift(FFTW.rfft(trace(t)), 2)
    delta_freq = 1/(nsamples(t)*t.delta)
    FourierDASArray{T,Matrix{Complex{TV}},P}(
        starttime(t), delta_freq, t.starting_distance, t.distance_spacing,
        nsamples(t), t.sta, t.evt, data, t.picks, t.meta
    )
end

function Seis.ifft(
    f::FourierDASArray{T,<:AbstractMatrix{<:Complex{TM}},P},
    d::Union{Integer,Nothing}=nothing
) where {T,TM,P}
    # Use original number of points if the length of the frequency coefficients
    # array has not been changed
    n = if f.nsamples in 2*nfrequencies(f) .- (1, 2)
        f.nsamples
    # Otherwise, we can't know how long the time-domain trace is meant
    # to be and so assume we want an even number of points unless we are told
    # how many we should have
    elseif d === nothing
        2*nfrequencies(f) - 2
    else
        if d == 2*nfrequencies(f) - 1 || d == 2*nfrequencies(f) - 2
            d
        else
            throw(ArgumentError("d must be either 2nf - 2 or 2nf - 1, " *
                "where nf is the number of frequency points"))
        end
    end

    # Scale the amplitude appropriately in case we have changed the
    # spectrum length
    scale = n/f.nsamples
    # Ensure we pass `d` as an `Int` (not `Int64`) on x86, otherwise there
    # is no method to construct a `FFTW.FakeArray`.
    data = FFTW.irfft(FFTW.ifftshift(trace(f), 2), Int(n)).*scale

    if size(data, 2) != length(f.sta)
        error(
            "number of stations ($(length(f.sta))) does not agree with " *
            "size of inverse-transformed matrix ($(size(data, 2))).  " *
            "Changing channel count is currently unsupported in the " *
            "Fourier domain."
        )
    end

    delta = 1/(n*f.delta)

    DASArray{T,Matrix{TM},P}(
        Seis.starttime(f), delta, f.starting_distance, f.distance_spacing,
        f.evt, f.sta, data, f.picks, f.meta
    )
end

"""
    nwavenumbers(f::FourierDASArray) -> n

Return the number of wavenumbers in the frequency-domain DAS array `f`
which are distinct only in terms of magnitude.  For a trace array with
original length `N`, this is `N÷2 + 1`.
"""
nwavenumbers(f::FourierDASArray) = size(trace(f), 2)÷2 + 1

"""
    wavenumbers(f::FourierDASArray)

Return the wavenumbers (spatial frequencies) of the Fourier-domain
DAS array `f` for each point of the second dimension of the underlying
trace data obtained with `trace(f)`.  The units are m⁻¹, i.e., there
is no factor of 2π.
"""
function wavenumbers(f::FourierDASArray)
    n = length(f.sta)
    if iseven(n)
        ((-n÷2):(n÷2 - 1))./(n*f.distance_spacing)
    else
        ((-(n - 1)÷2):((n - 1)÷2))./(n*f.distance_spacing)
    end
end
