using Seis
using TraceArrays
using Test

@testset "Operations" begin
    @testset "integrate" begin
        ntraces = 4
        ts = [Trace(1.5, 0.1, rand(20)) for _ in 1:ntraces]

        @testset "$T" for T in (TraceArray, DASArray)
            arr = if T == TraceArray
                T(ts)
            else
                T(ts, 5, 10)
            end

            @testset "Default" begin
                @test integrate(arr) == integrate(arr, :trapezium)
            end

            @testset "In-place" begin
                @test integrate(arr, :trapezium) == integrate!(deepcopy(arr), :trapezium)
                @test integrate(arr, :rectangle) == integrate!(deepcopy(arr), :rectangle)
            end

            @testset "Errors" begin
                @testset "Method" begin
                    @test_throws ArgumentError integrate(arr, :unsupported_method)
                    @test_throws ArgumentError integrate!(arr, :unsupported_method)
                end

                @testset "Length" begin
                    @test_throws ArgumentError integrate(cut(arr, times(arr)[[1,1]]...), :trapezium)
                end
            end

            @testset "Single-trace test" begin
                @testset "$method" for method in (:trapezium, :rectangle)
                    for i in 1:ntraces
                        @test integrate(arr[i], method) == integrate(arr, method)[i]
                    end
                end                
            end
        end
    end
end
