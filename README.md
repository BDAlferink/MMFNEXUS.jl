# MMFNEXUS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BDAlferink.github.io/MMFNEXUS.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BDAlferink.github.io/MMFNEXUS.jl/dev/)
[![Build Status](https://github.com/BDAlferink/MMFNEXUS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BDAlferink/MMFNEXUS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/BDAlferink/MMFNEXUS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BDAlferink/MMFNEXUS.jl)


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

The basic usage of the NEXUS+ routine requries a density field with non-zero values everywhere. For optimal results, we suggest density field reconstructions using [DTFE](https://github.com/MariusCautun/DTFE), or [Phase-Space DTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl). How to reconstruct a density field from a particle distribution is found on the respective pages. The density field should be normalized, i.e. $\frac{\rho}{\rho_{\text{mean}}} = 1 + \delta$.

Give the normalized density field (`densityField`), we identify the cosmic web environments as follows:

```julia
using MMFNEXUS

# field and box parameters (example from the reconstruction of the illustris-3 box sampled at 256^3)
N = 256 # number of gridpoints per dimension
L = 75. # Box size in cMpc/h
totalMass = 4e8 * 455^3 # total mass contained in simulation box in Msun/h 

MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(densityField, N, L, totalMass);

```

The resulting `MMF_*` outputs are BitArray's of size (N,N,N) where for each voxel, one of the morphological environments has the value `1` and all others are `0` to indicate to which environment it belongs. There are a number of optinal settings which are discussed in the tutorial section of the documentation. A theory section is also found there.

## Contributors
This Julia implementation is written by:
- Bram Alferink ([alferink@astro.rug.nl](mailto:alferink@astro.rug.nl))

The original NEXUS+ algorithm published as [NEXUS: tracing the cosmic web connection (Marius Cautun , Rien van de Weygaert , Bernard J. T. Jones)](https://academic.oup.com/mnras/article/429/2/1286/1038906) is written by:
- Marius Cautun

We thank:
- Ivan Spirov
- Rien van de Weijgaert
- Job Feldbrugge