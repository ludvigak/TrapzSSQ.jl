using Test
using TrapzSSQ
using QuadGK

@testset "TrapzSSQ" begin
    include("test_pow.jl")
    include("test_log.jl")
    include("test_weights.jl")
    include("test_interface.jl")
    include("test_starfish.jl")
end
;
