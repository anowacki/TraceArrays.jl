using Dates: now
using Seis
using Test
using TraceArrays

import FFTW

@testset "FourierDASArray" begin
    @testset "Construction errors" begin
        @test_throws ArgumentError FourierDASArray(
            b=0, delta=1, nsamples=-1, starting_distance=0, distance_spacing=1,
            data=rand(ComplexF64, 10, 10)
        )
        @test_throws ArgumentError FourierDASArray(
            b=0, delta=Inf, nsamples=21, starting_distance=0, distance_spacing=1,
            data=rand(ComplexF64, 10, 10)
        )
    end

    @testset "Type $T" for T in (Float16, Float32, Float64, BigFloat)
        @testset "Matrix eltype $TM" for TM in (Float32, Float64)
            M = Matrix{Complex{TM}}
            @testset "Position $P" for P in (Seis.Cartesian{T}, Seis.Geographic{T})
                @testset "Construction" begin
                    b = 100*randn(T)
                    delta = rand(T)/10 # Time-domain sampling interval
                    starting_distance = 10*randn(T)
                    distance_spacing = 10*rand(T)
                    picks = Seis.SeisDict{Union{Int,Symbol},Seis.Pick{T}}(
                        Dict(:A=>Seis.Pick{T}(1.0, "A"))
                    )
                    meta = Seis.SeisDict{Symbol,Any}(
                        Dict(:label=>"text")
                    )
                    npts = rand(20:41)
                    nchan = rand(20:41)
                    evt = Event{T,P}(time=now())
                    sta = [Station{T,P}(sta=string(i)) for i in 1:nchan]
                    data = randn(Complex{TM}, npts, nchan)
                    
                    @testset "Default" begin
                        t = FourierDASArray{T,M,P}(;
                            b, delta, starting_distance, distance_spacing,
                            evt, sta, picks, meta, data, nsamples=npts,
                        )
                        @test t isa FourierDASArray{T,M,P}
                        @test t.b == b
                        @test t.delta == 1/(npts*delta)
                        @test t.starting_distance == starting_distance
                        @test t.distance_spacing == distance_spacing
                        @test t.meta == meta
                        @test trace(t) == data
                        @test t == FourierDASArray{T,M,P}(
                            b, 1/(npts*delta), starting_distance, distance_spacing,
                            npts, sta, evt, data, picks, meta
                        )
                    end
                end

                @testset "fft" begin
                    t = DASArray(;
                        b=rand(T), delta=rand(T),
                        starting_distance=-20, distance_spacing=0.5,
                        data=rand(TM, 40, 20),
                        T, M=Matrix{TM}, P
                    )
                    f = fft(t)

                    @testset "Fields carried over" begin
                        @testset "Field $field" for field in (
                            :b, :starting_distance, :distance_spacing,
                            :sta, :evt, :picks, :meta
                        )
                            @test getfield(t, field) == getfield(f, field)
                        end
                    end

                    @testset "Fields computed" begin
                        @test f.nsamples == nsamples(t)
                        @test f.delta == 1/(f.nsamples*t.delta)
                    end
                end

                @testset "fft-ifft round-trip" begin
                    @testset "nsamples $odd_even" for odd_even in (:odd, :even)
                        npts = 40 + (odd_even == :odd)
                        t = DASArray(;
                            b=rand(T), delta=rand(T),
                            starting_distance=-20, distance_spacing=0.5,
                            data=rand(TM, npts, 20),
                            T, M=Matrix{TM}, P
                        )
                        t′ = ifft(fft(t))
                        @test length(t) == length(t′)
                        @test length(t′.sta) == size(trace(t′), 2)
                        @testset "Field $f" for f in fieldnames(typeof(t))
                            if f === :data
                                @test t.data ≈ t′.data
                            elseif f === :delta
                                @test t.delta ≈ t′.delta rtol=sqrt(eps(T))
                            else
                                @test getfield(t, f) == getfield(t′, f)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Getters" begin
        # Ensure we check both odd and even numbers of samples and channels
        @testset "nsamples/nchannels + $plus" for plus in (0, 1)
            npts, nchan = rand(20:41) + plus, rand(20:41) + plus
            delta = rand() + eps()
            distance_spacing = rand()
            t = DASArray(;
                b=rand(), delta, starting_distance=rand(),
                distance_spacing, data=randn(npts, nchan)
            )
            f = fft(t)
            @test collect(frequencies(f)) ≈ FFTW.fftshift(FFTW.rfftfreq(npts, 1/delta))
            @test nfrequencies(f) == length(FFTW.rfftfreq(npts, 1/delta))
            @test collect(wavenumbers(f)) ≈ FFTW.fftshift(FFTW.fftfreq(nchan, 1/distance_spacing))
        end
    end

    @testset "Synthetic case" begin
        v_app = 5 # Apparent velocity m/s;
        freq = 0.5 # Wave frequency, Hz;
        delta = 0.05 # Sampling interval;
        distance_spacing = 0.05 # Channel spacing, m;
        # Sine wave of frequency 1 Hz with apparent velocity of 5 m/s
        t = DASArray(;
            b=0, delta, starting_distance=0, distance_spacing,
            data=sin.(2π*freq.*((0:201).*delta .- (0:201)'.*distance_spacing./v_app))
        )
        T = fft(t);

        # Location of greatest magnitude
        indmax = argmax(abs.(trace(T)))
        @test indmax == CartesianIndex(6, 101)
        # Peak frequency as input
        f = frequencies(T)[indmax[1]]
        @test f ≈ freq rtol=0.01
        # Negative wavenumbers are for positive-travelling waves
        λ⁻¹ = -wavenumbers(T)[indmax[2]]
        # Apparent velocity is as we put in
        @test f/λ⁻¹ ≈ v_app rtol=0.01
    end

    @testset "Constant channel count requirement" begin
        t = DASArray(b=0, delta=1, starting_distance=-100,
            distance_spacing=0.1, data=rand(Float32, 20, 20)
        )
        f = fft(t)
        f.data = trace(f)[:,2:end-1]
        @test_throws ErrorException ifft(f)
        f.data = FFTW.rfft(trace(t))
        @test ifft(f) isa DASArray
    end
end
