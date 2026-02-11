using HDF5
using Statistics:mean
using MMFNEXUS
using JLD2

function load_data(file)
    fid = h5open(file, "r")
    densityfield = read(fid["densityfield"])
    close(fid)

    densityfield = densityfield ./ mean(densityfield)
    return densityfield
end
densityfield = load_data("test/data/densityfield.h5");


filter_scales = 6 #(0:1)
# filter_scales = (0.5, .5^.5)

density_contrast_node = 370.
min_node_mass = 1e13 #Msun/h
min_fila_volume = 10 #(Mpc/h)^3
min_wall_volume = 10 #(Mpc/h)^3
min_scale = .5 #minimum smoothing scale in Mpc/h

N = 64 # number of gridpoints per dimension
L = 50. # Box size in cMpc/h
M = 4.075161606358443e10 * 64^3 # total mass contained in the box in in Msun


# using MMFNEXUS
@time MMF_node, MMF_fila, MMF_wall, MMF_void = NEXUS_Plus(densityfield, N, 
L, M, filter_scales, density_contrast_node, min_node_mass, min_fila_volume, min_wall_volume; R0 = min_scale, level = :debug);

test_output = NEXUS_Plus(densityfield, N, 
L, M, filter_scales, density_contrast_node, min_node_mass, min_fila_volume, min_wall_volume; R0 = min_scale, level = :debug);
@save "test/data/NEXUSTestClassification.jld2" test_output
typeof(test_output)
println("Done with MMF-NEXUS classification!")

