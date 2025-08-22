using TraceArrays
using Test

@testset "TraceArrays.jl" begin
    include("tracearray.jl")
    include("dasarray.jl")
    include("fourierdasarray.jl")
    include("operations.jl")
end
