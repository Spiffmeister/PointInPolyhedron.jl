using Test
# using SPECreader

@testset "Unit cube" begin
    include("unitcube.jl")
end

@testset "Torus" begin
    include("torus.jl")
end
