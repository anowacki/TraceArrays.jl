using Seis
using Test
using TraceArrays

@testset "TraceArray" begin
    @testset "append!" begin
        @testset "Errors" begin
            @testset "Wrong nsamples" begin
                @test_throws ArgumentError TraceArray(sample_data(:array))
            end

            @testset "Empty" begin
                @test_throws ArgumentError TraceArray(typeof(sample_data())[])
            end
        end
    end

    @testset "AbstractArray interface" begin
        @testset "Base.eltype" begin
            # FIXME: eltype of a TraceArray should be `Trace` after changes to
            #        getindex.
            @test eltype(TraceArray(b=10, delta=0.1, data=rand(Float16, 2, 2))) == Float16
        end
        @testset "Base.pairs" begin
            for (i, t) in pairs(TraceArray(b=0, delta=1, data=Float64[i for j in 1:2, i in 1:4]))
                @test t == Seis.Trace(b=0, delta=1, data=fill(i, 2))
            end
        end

        @testset "Base.keys" begin
            @test keys(TraceArray(b=0, delta=1, data=rand(3, 4))) == 1:4
        end

        @testset "Base.iterate" begin
            t = TraceArray(b=0, delta=1, data=Float64[1 2 3; 2 3 4])
            for (i, tt) in pairs(t)
                @test tt == Seis.Trace(b=0, delta=1, data=Seis.trace(t)[:,i])
            end
        end

        @testset "Base.collect" begin
            t = TraceArray(b=rand(), delta=rand(), data=rand(4, 3))
            @test collect(t) == [tt for tt in t]
        end
    end
end
