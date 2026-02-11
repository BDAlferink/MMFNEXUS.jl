using MMFNEXUS
using Test
using HDF5
using Statistics:mean
using JLD2

# Load test density field
function load_data(file)
    fid = h5open(file, "r")
    densityfield = read(fid["densityfield"])
    close(fid)

    densityfield = densityfield ./ mean(densityfield)
    return densityfield
end

@testset "MMFNEXUS.jl" begin

    densityfield = load_data("data/densityfield.h5");

    # Load test comparison data
    @load "data/NEXUSTestClassification.jld2" test_output

    # set NEXUS+ parameters
    min_scale = .5 #minimum smoothing scale in Mpc/h
    filter_scales = 6 #max n in min_scale*(√2)^n, starting at n=0
    density_contrast_node = 370.
    min_node_mass = 1e13 #Msun/h
    min_fila_volume = 10 #(Mpc/h)^3
    min_wall_volume = 10 #(Mpc/h)^3

    # set box parameters
    N = 64 # number of gridpoints per dimension
    L = 50. # Box size in Mpc/h
    M = 4.075161606358443e10 * 64^3 # total mass contained in the box in in Msun

    NEXUSTest = NEXUS_Plus(densityfield, N, L, M, filter_scales, density_contrast_node, min_node_mass, min_fila_volume, min_wall_volume; R0 = min_scale, level = :none);

    #test nodes
    @test NEXUSTest[1] ≈ test_output[1]
    #test filaments
    @test NEXUSTest[2] ≈ test_output[2]
    #test walls
    @test NEXUSTest[3] ≈ test_output[3]
    #test voids
    @test NEXUSTest[4] ≈ test_output[4]


end
