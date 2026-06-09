module MMFNEXUS

    # ---- imports ----
    using FFTW
    using Images
    using Plots
    using Statistics:mean
    using Roots
    using Base.Threads
    using Logging
    using ImageFiltering

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
