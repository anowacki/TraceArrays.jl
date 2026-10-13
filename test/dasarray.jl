using Test
using TraceArrays

@testset "DASArrays" begin
    @testset "Integration" begin
        @testset "_integrate N = $N" for N in 5:5:20
            @testset "n = $n" for n in 3:2:7
                if N < n
                    continue
                end
                spacing = rand(2:10)
                u = rand(N)
                U = similar(u, (N - n + 1))
                @test TraceArrays._integrate!(U, u, spacing/2, n) ≈
                    [spacing*(u[i]/2 + u[i+n-1]/2 + sum(u[(i+1):(i+n-2)])) for i in eachindex(u) if i <= N - n + 1]
            end
        end
    end

    @testset "Conversion to TraceArray" begin
        t = DASArray(;
            b=1, delta=0.1, starting_distance=10, distance_spacing=5, data=randn(4, 3)
        )
        t′ = TraceArray(t)

        @test t′ isa TraceArray
        @test t′.delta == t.delta
        @test t′.b == t.b
        @test t′.evt == t.evt
        @test t′.sta == t.sta
        @test t′.data == t.data
    end
end