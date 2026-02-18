```@meta
CurrentModule = MMFNEXUS
```

# MMF NEXUS

Documentation for [MMFNEXUS](https://github.com/BDAlferink/MMFNEXUS.jl).

UNDER DEVELOPMENT

This is the Julia implementation for the NEXUS+ algorithm. More information and standalone application can be found on our website in the future.

## Installation

<!--
The Multiscale Morphology Filter (MMF) NEXUS can be installed with the Julia package manager. From the Pkg Repl mode run:

add pkg thing in julia
-->

For now, the Multiscale Morphology Filter (MMF) NEXUS can not be installed through the package manager yet. This will become available soon.

## Usage

For now, the NEXUS+ algorithm is the only routine available. This requries a density field with non-zero values everywhere. For optimal results, we suggest density field reconstructions using [DTFE](https://github.com/MariusCautun/DTFE), or [Phase-Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl). On how to reconstruct a density field from a particle distribution is found on the respective pages. The density field should be normalized, i.e. $\frac{\rho}{\rho_{\text{mean}}} = 1 + \delta$.

Give the normalized density field (`densityField`), we identify the cosmic web environments as follows:

```julia
using MMFNEXUS

# field and box parameters (example from the reconstruction of the illustris-3 box sampled at 256^3)
N = 256 # number of gridpoints per dimension
L = 75. # Box size in cMpc/h
totalMass = 4e8 * 455^3 # total mass contained in simulation box in Msun/h 

# NEXUS+ parameters
minimumFilterScale = 0.5 # minimum smoothing scale in Mpc/h
filter_scales = 6 # maximum index n in min_scale*(√2)^n. Here it does n = 0,1,2,3,4,5,6
density_contrast_node = 370.
minimum_node_mass = 1e13 # Msun/h
minimum_filament_volume = 5 # (MPC/h)^3
minimum_wall_volume = 5 # (MPC/h)^3

MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(densityField, N, L, totalMass, filter_scales, density_contrast_node, minimum_node_mass, minimum_filament_volume, minimum_wall_volume; R0 = min_scale);

```

The resulting `MMF_*` outputs are BitArray's of size (N,N,N) where for each voxel, one of the morphological environments has the value `1` and all others are `0` to indicate to which environment it belongs.

For more details, theory, and a tutorial, please consult the Documentation.

## Contributors
This Julia implementation is written by:
- Bram Alferink ([alferink@astro.rug.nl](mailto:alferink@astro.rug.nl))

The original NEXUS+ algorithm published as [NEXUS: tracing the cosmic web connection (Marius Cautun , Rien van de Weygaert , Bernard J. T. Jones)](https://academic.oup.com/mnras/article/429/2/1286/1038906) is written by:
- Marius Cautun

We thank:
- Ivan Spirov
- Rien van de Weijgaert
- Job Feldbrugge