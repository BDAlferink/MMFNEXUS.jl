```@meta
CurrentModule = MMFNEXUS
```

# Tutorial
In this tutorial, we demonstrate how to use MMF NEXUS on a density field. The density field used in this example is obtained from using [Phase Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl) on the particle distributino of a 64^3 GADGET-4 simulation. 

In principle MMF NEXUS can be applied to any continous density field. Due to the filtering in log space, NEXUS+ requires the field to be positive valued everywhere. We suggest to use [DTFE](https://github.com/MariusCautun/DTFE) or [Phase Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl) for the density field reconstructions as those methods preserve the geometric properties of the matter distribution and are positive valued at each point.

We now start by importing relevant libraries and loading the data:

```@example tutorial1
using MMFNEXUS, HDF5, Statistics

# Load test density field
function load_data(file)
    fid = h5open(file, "r")
    densityfield = read(fid["densityfield"])
    close(fid)

    densityfield = densityfield ./ mean(densityfield)
    return densityfield
end

densityfield = load_data("assets/data/densityfield.h5");
nothing

```
lets do a plot here

```@example tutorial1
# set up the box
N = 64
L = 50. #Mpc/h
M = 4.075161606358443e10 * 64^3 # total mass contained in the box in in Msun

# set NEXUS+ parameters
min_scale = .5 #minimum smoothing scale in Mpc/h
filter_scales = 6 #max n in min_scale*(√2)^n, starting at n=0
density_contrast_node = 370.
min_node_mass = 1e13 #Msun/h
min_fila_volume = 10 #(Mpc/h)^3
min_wall_volume = 10 #(Mpc/h)^3
nothing

```
calculation here.
etc.

## Options


## Multithreading


