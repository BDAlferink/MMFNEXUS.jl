```@meta
CurrentModule = MMFNEXUS
```

# Tutorial
In this tutorial, we demonstrate how to use MMF NEXUS on a density field. The density field used in this example is obtained from using [Phase Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl) on the particle distributino of a 64^3 GADGET-4 simulation. 

In principle MMF NEXUS can be applied to any continous density field. Due to the filtering in log space, NEXUS+ requires the field to be positive valued everywhere. We suggest to use [DTFE](https://github.com/MariusCautun/DTFE) or [Phase Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl) for the density field reconstructions as those methods preserve the geometric properties of the matter distribution and are positive valued at each point.

We now start by importing relevant libraries and loading the data. We plot a slice of the density field to illustrate what we are working with.

```@example tutorial1
using MMFNEXUS, HDF5, Statistics, Plots

# Load test density field
function load_data(file)
    fid = h5open(file, "r")
    densityfield = read(fid["densityfield"])
    close(fid)

    densityfield = densityfield ./ mean(densityfield)
    return densityfield
end

densityfield = load_data("assets/data/densityfield.h5");
heatmap(log10.(densityfield[:, :, 20]), aspect_ratio=:equal, c=:viridis, title="Density Field Slice", xlabel="x [Mpc/h]", ylabel="y [Mpc/h]", colorbar_title="log₁₀(1+δ)", xlims=(0, 50), ylims=(0, 50), xticks=0:10:50, yticks=0:10:50);


```
We now set up the simulation box.

```@example tutorial1
# set up the box
N = 64 # number of voxels per side
L = 50. # side length in Mpc/h
M = 4.075e10 * 64^3 # total mass contained in the box in in Msun

```
The NEXUS+ routine is called with the function `NEXUS_Plus`, which takes the density field and a number of keyword arguments. These will be explained below. One of the keyword arguments is a verbose level which can be set to `none`, `info`, or `debug`. `info` gives some basic information on what the calulation is doing. `debug` gives additional information and figures of the threshold calulation to see what is going on. `NEXUS_Plus` gives four boolean matrices where for each environemt a `1` means the voxel is considered to be in the corresponding environemnt

```@example tutorial1
MMF_node, MMF_fila, MMF_wall, MMF_void = NEXUS_Plus(densityfield, N, L, M; filter_parse = 6 ,R0 = .5, level = :info);
nothing
```

We can use the result to plot the contours on top of the slice of a density field, which we show below, or any subsequent analysis.

```@example tutorial1
slice_index = 20
heatmap(log10.(densityfield[:, :, slice_index]), aspect_ratio=:equal, c=:viridis, title="Density Field Slice", xlabel="x [Mpc/h]", ylabel="y [Mpc/h]", colorbar_title="log₁₀(1+δ)", xlims=(0, 50), ylims=(0, 50), xticks=0:10:50, yticks=0:10:50)
contour!(MMF_wall[:, :, slice_index], levels=[0.5], color=:green, linewidth=2, label="Walls")
contour!(MMF_fila[:, :, slice_index], levels=[0.5], color=:blue, linewidth=2, label="Filaments")
contour!(MMF_node[:, :, slice_index], levels=[0.5], color=:red, linewidth=2, label="Nodes")

```

## Options

A number of parameters can be changed manually in NEXUS+. Particularly the minimum filtering scale can be changed according to the resolution of the simulation and corresponding density field. The keyword arguments and their defaults of ```NEXUS_Plus``` are: 

```julia
Δ::Real = 370. # density contrast for node detection
min_node_mass::Real = 1e13 # minimum mass of a node in Msun/h
min_fila_volume::Real = 10 # minimum volume of a filament in (Mpc/h)^3
min_wall_volume::Real = 10 # minimum volume of a wall in (Mpc/h)^3
R0::Real = 1. #minimum smoothing scale in Mpc/h
filter_parse = 4 #max n in min_scale*(√2)^n, starting at n=0
level::Symbol = :info # verbose level
```
Δ is the density contrast for node detection, this is further explained in the theory section. The minimum node mass is used as a physical selection criteria to avoid spurious detections of tiny dense clumps as clusters. Similarly,the minimum fila/wall argument is used to avoid spurious detections. R0 is the minimum filtering scale. Filter parse takes either an integer, in which case the filter scales are R0(√2)^n where n are integers from 0 up to the given integer, or a tuple/array of all the user defined scales in Mpc/h can be given directly. In case a tuple is given, the argument R0 does not do anything. The level parameter lets us select the verbose level which can be set to `none`, `info`, or `debug`. `info` gives some basic information on what the calulation is doing. `debug` gives additional information and figures of the threshold calulation to see what is going on.

The function with all keyword arguments can be called as:
```julia
MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(densityField, N, L, totalMass; filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, minimum_wall_volume = min_wall_volume; R0 = R0, level = level);
```

## Multithreading

Multithreading is `automatically` enabled. It will utilize the cores available to your instance of Julia. This can be set by opening julia with a specifified number of threads:

```bash
$ julia --threads 4
```
to open julia with 4 threads for example or by setting an environment variable. For more information, see the corresponding [multi-threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/).
