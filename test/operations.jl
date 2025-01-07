using Seis
using TraceArrays
using Test

@testset "Operations" begin
    @testset "differentiate" begin
        ntraces = 3
        ts = [Trace(0.34, 0.25, rand(20)) for _ in 1:ntraces]

        @testset "$T" for T in (TraceArray, DASArray)
            arr = if T == TraceArray
                T(ts)
            else
                T(ts, 7.5, 1.5)
            end

            @testset "Default" begin
                @test differentiate(arr) == differentiate(arr; points=2)
            end

            @testset "In-place" begin
                @test differentiate(arr) == differentiate!(deepcopy(arr))
            end

            @testset "Errors" begin
                @testset "points" begin
                    @testset "$func" for func in (differentiate, differentiate!)
                        @test_throws ArgumentError func(arr; points=0)
                        @test_throws ArgumentError func(arr; points=1)
                        @test_throws ArgumentError func(arr; points=4)
                        @test_throws ArgumentError func(arr; points=6)
                    end
                end

                @testset "Length of trace" begin
                    @testset "$points points" for points in (2, 3, 5)
                        @testset "$func" for func in (differentiate, differentiate!)
                            short_arr = TraceArray(
                                [Trace(-0.3, 0.1, randn(points - 1)) for _ in 1:ntraces]
                            )
                            @test_throws ArgumentError func(short_arr; points)
                        end
                    end
                end
            end

            @testset "Single-trace test" begin
                @testset "$points" for points in (2, 3, 5)
                    for i in 1:ntraces
                        @test differentiate(arr[i]; points) == differentiate(arr; points)[i]
                    end
                end                
            end

        end
    end

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
