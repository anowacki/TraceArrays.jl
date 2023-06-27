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
end