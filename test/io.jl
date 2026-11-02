using Dates: DateTime
import Seis
using Test
using TraceArrays

@testset "IO" begin
    test_data_path(file...) = joinpath(@__DIR__, "data", file...)

    @testset "OptoDAS" begin
        file = test_data_path("optodas.hdf5")
        @testset "Plain reading" begin
            t = read_optodas(file)
            @test t isa DASArray
            @test t.delta == 0.01
            @test distances(t) == [-1, 9]
            @test t.meta.gauge_length == 20.42600381549363
            @test Seis.origin_time(t) == DateTime("2023-02-06T11:12:41")
            @test Seis.trace(t) == [-1 10; -2 20; -3 30; -4 40]
        end

        @testset "xdecimate" begin
            @test distances(read_optodas(file; xdecimate=2)) == [-1]
        end

        @testset "Distance range" begin
            @test distances(read_optodas(file; xlim=(-2, -1))) == [-1]
            @test distances(read_optodas(file; xlim=(-0.99, Inf))) == [9]
        end

        @testset "Time range" begin
            @testset "Times" begin
                t1 = read_optodas(file; tlim=(-Inf, 0.02))
                @test Seis.trace(t1) == [-1 10; -2 20; -3 30]
                @test starttime(t1) == 0.0

                t2 = read_optodas(file; tlim=(0.001, 0.02))
                @test Seis.trace(t2) == [-2 20; -3 30]
                @test starttime(t2) == 0.01
            end

            @testset "Dates" begin
                t1 = read_optodas(file; tlim=(DateTime(1970), DateTime("2023-02-06T11:12:41.015")))
                @test starttime(t1) == 0.00
                @test Seis.trace(t1) == [-1 10; -2 20]

                t1 = read_optodas(file; tlim=(DateTime("2023-02-06T11:12:41.005"), DateTime("2023-02-06T11:12:41.03")))
                @test starttime(t1) == 0.01
                @test Seis.trace(t1) == [-2 20; -3 30; -4 40]
            end
        end

        @testset "Header only" begin
            t_full = read_optodas(file)
            t_header = read_optodas(file; header_only=true)
            @test typeof(t_full) == typeof(t_header)

            @testset "Field $field" for field in fieldnames(typeof(t_full))
                value_full = getfield(t_full, field)
                value_header = getfield(t_header, field)
                if field == :data
                    @test isempty(value_header)
                    @test size(value_header) == (0, 0)
                else
                    @test value_full == value_header
                end
            end
        end

        @testset "_epoch_to_datetime" begin
            @test TraceArrays.InputOutput._epochseconds_to_datetime(0) == DateTime(1970)
            @test TraceArrays.InputOutput._epochseconds_to_datetime(3661.123) == DateTime(1970, 1, 1, 1, 1, 1, 123)
            @test(
                @test_logs(
                    (
                        :warn,
                        "Unexpected non-integer number of milliseconds since Unix epoch in OptoDAS header"
                    ),
                    TraceArrays.InputOutput._epochseconds_to_datetime(3661.1235)
                ) == DateTime(1970, 1, 1, 1, 1, 1, 124)
            )
        end

        @testset "_cut_indices" begin
            @test TraceArrays.InputOutput._cut_indices(-80:10:20, (-79, 0)) == 2:9
            @test TraceArrays.InputOutput._cut_indices(-80:10:20, (-70, 1)) == 2:9
        end
    end
end
