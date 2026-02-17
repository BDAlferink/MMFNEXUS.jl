using HDF5
using Statistics:mean
using MMFNEXUS
using JLD2
using Plots



function load_data(file)
    fid = h5open(file, "r")
    densityfield = read(fid["densityfield"])
    close(fid)

    densityfield = densityfield ./ mean(densityfield)
    return densityfield
end
densityfield = load_data("data/densityfield.h5");

# Plot a slice of the density field
heatmap(log10.(densityfield[:, :, 20]), 
    aspect_ratio=:equal, 
    c=:viridis, 
    title="Density Field Slice",
    xlabel="x [Mpc/h]", 
    ylabel="y [Mpc/h]",
    colorbar_title="log₁₀(1+δ)",
    xlims=(0, 50),
    ylims=(0, 50),
    xticks=0:10:50,
    yticks=0:10:50)
# savefig("densityfield_slice.png")

filter_scales = 4 #(0:1)
# filter_scales = (0.5, .5^.5)

density_contrast_node = 370.
min_node_mass = 1e13 #Msun/h
min_fila_volume = 10 #(Mpc/h)^3
min_wall_volume = 10 #(Mpc/h)^3
min_scale = 1. #minimum smoothing scale in Mpc/h

N = 64 # number of gridpoints per dimension
L = 50. # Box size in cMpc/h
M = 4.075e10 * 64^3 # total mass contained in the box in in Msun


# using MMFNEXUS
@time MMF_node, MMF_fila, MMF_wall, MMF_void = NEXUS_Plus(densityfield, N, 
L, M, filter_scales, density_contrast_node, min_node_mass, min_fila_volume, min_wall_volume; R0 = min_scale, level = :info);

function plot_cosmic_web_slice(densityfield, MMF_wall, MMF_fila, MMF_node, idx)
    heatmap(log10.(densityfield[:, :, idx]), 
        aspect_ratio=:equal, 
        c=:viridis, 
        title="Density Field with Cosmic Web",
        xlabel="x [Mpc/h]", 
        ylabel="y [Mpc/h]",
        colorbar_title="log₁₀(1+δ)",
        xlims=(0, 50),
        ylims=(0, 50),
        xticks=0:10:50,
        yticks=0:10:50)
    
    contour!(MMF_wall[:, :, idx], levels=[0.5], color=:green, linewidth=2, label="Walls")
    contour!(MMF_fila[:, :, idx], levels=[0.5], color=:blue, linewidth=2, label="Filaments")
    contour!(MMF_node[:, :, idx], levels=[0.5], color=:red, linewidth=2, label="Nodes")
end

plot_cosmic_web_slice(densityfield, MMF_wall, MMF_fila, MMF_node, 20)

test_output = NEXUS_Plus(densityfield, N, 
L, M, filter_scales, density_contrast_node, min_node_mass, min_fila_volume, min_wall_volume; R0 = min_scale, level = :debug);
@save "test/data/NEXUSTestClassification.jld2" test_output
typeof(test_output)
println("Done with MMF-NEXUS classification!")

