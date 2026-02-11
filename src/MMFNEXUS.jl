module MMFNEXUS

    # ---- imports ----
    using HDF5
    using FFTW
    using LinearAlgebra
    using Images
    using Statistics:mean
    using Roots
    using Base.Threads
    using Logging

    # ---- set threads ----
    FFTW.set_num_threads(Threads.nthreads())

    # ---- includes ----
    # NEXUS+ routine
    include("NEXUSPlus.jl")

    # utilities for logging, input and math
    include("utils.jl")

    # signature calculations
    include("response.jl")

    # threshold calculations
    include("threshold.jl")

    export NEXUS_Plus

    
end
