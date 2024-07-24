using TraceArrays
using Test

@testset "TraceArrays.jl" begin
    include("tracearray.jl")
    include("dasarray.jl")
    include("operations.jl")
end
