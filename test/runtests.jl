using MMFNexus
using Test
using HDF5
using Statistics:mean
using JLD2

# Load test density field
function load_data(file)
    fid = h5open(file, "r") do fid
        densityfield = read(fid["densityfield"])
        return densityfield ./ mean(densityfield)
    end
end

const densityfield = load_data("data/densityfield.h5")
const N = 64 # number of gridpoints per dimension
const L = 50. # Box size in cMpc/h
const M = 4.075161606358443e10 * 64^3 # total mass contained in the box in Msun

const KWARGS = (
    filter_parse = 6,
    Δ = 370.0,
    min_node_mass = 1e13,
    min_fila_volume = 10,
    min_wall_volume = 10,
    R0 = 0.5,
    level = :none,
)

@testset "MMFNexus.jl" begin

    @testset "SimBox" begin
        sim = MMFNexus.SimBox(N, L, M)
        @test sim.ρ_mean ≈ M / L^3
        @test sim.voxel_volume ≈ (L / N)^3
        @test sim.N == N
    end

    @testset "parse_filter_scales" begin
        # integer dispatch
        scales = MMFNexus.parse_filter_scales(4, 0.5)
        @test length(scales) == 5          # n = 0:4
        @test scales[1] ≈ 0.5
        @test scales[2] ≈ 0.5 * sqrt(2)

        # vector dispatch passthrough
        explicit = [0.5, 1.0, 2.0]
        @test MMFNexus.parse_filter_scales(explicit) == explicit
    end

    @testset "compute_eigenvalues_sym3" begin
        # diagonal matrix — known eigenvalues
        l1, l2, l3 = MMFNexus.compute_eigenvalues_sym3(3.0, 1.0, 2.0, 0.0, 0.0, 0.0)
        @test l1 ≈ 1.0
        @test l2 ≈ 2.0
        @test l3 ≈ 3.0

        # known symmetric matrix with known eigenvalues
        # [ 2  1  0 ]  eigenvalues: 1, 1, 3
        # [ 1  2  0 ]
        # [ 0  0  1 ]
        l1, l2, l3 = MMFNexus.compute_eigenvalues_sym3(2.0, 2.0, 1.0, 1.0, 0.0, 0.0)
        @test l1 ≈ 1.0 atol=1e-10
        @test l2 ≈ 1.0 atol=1e-10
        @test l3 ≈ 3.0 atol=1e-10

        # output is always sorted ascending
        l1, l2, l3 = MMFNexus.compute_eigenvalues_sym3(1.0, 3.0, 2.0, 0.5, 0.3, 0.1)
        @test l1 ≤ l2 ≤ l3
    end    

    @testset "NEXUS_Plus method: Fourier" begin
        @load "data/NEXUSTestClassification.jld2" test_output
        result = NEXUS_Plus(densityfield, N, L, M; KWARGS..., method = :fourier)
        @test result[1] == test_output[1]   # nodes   — BitArray, use == not ≈
        @test result[2] == test_output[2]   # filaments
        @test result[3] == test_output[3]   # walls
        @test result[4] == test_output[4]   # voids
    end

    @testset "NEXUS_Plus: finitediff method" begin
        @load "data/NEXUSTestClassification_finitediff.jld2" test_output_finitediff
        result = NEXUS_Plus(densityfield, N, L, M; KWARGS..., method = :finitediff)
        @test result[1] == test_output_finitediff[1]
        @test result[2] == test_output_finitediff[2]
        @test result[3] == test_output_finitediff[3]
        @test result[4] == test_output_finitediff[4]
    end

    @testset "Invalid inputs" begin
        @test_throws ArgumentError NEXUS_Plus(densityfield, N+1, L, M; KWARGS...)           # wrong N
        @test_throws ArgumentError NEXUS_Plus(zeros(N,N,N), N, L, M; KWARGS...)             # non-positive ρ
        @test_throws ArgumentError NEXUS_Plus(densityfield, N, -1.0, M; KWARGS...)          # L ≤ 0
        @test_throws ArgumentError NEXUS_Plus(densityfield, N, L, -1.0; KWARGS...)          # M ≤ 0
        @test_throws Exception     NEXUS_Plus(densityfield, N, L, M; KWARGS..., method = :invalid)
        @test_throws Exception     NEXUS_Plus(densityfield, N, L, M; KWARGS..., level = :invalid)
    end

end
