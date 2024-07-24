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
end
