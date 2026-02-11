struct LevelMessageLogger <: AbstractLogger
    min_level::LogLevel
end

Logging.min_enabled_level(logger::LevelMessageLogger) = logger.min_level
Logging.shouldlog(logger::LevelMessageLogger, level, _module, group, id) =
    level ≥ logger.min_level

function Logging.handle_message(logger::LevelMessageLogger, level, message,
                                _module, group, id, file, line; kwargs...)
    println("[$(uppercase(string(level)))] $message")
end
function setup_logging(level::Symbol)
    logger = level === :none  ? NullLogger() :
             level === :info  ? LevelMessageLogger(Logging.Info) :
             level === :debug ? LevelMessageLogger(Logging.Debug) :
             error("Unknown log level")

    global_logger(logger)
end

function debugplot(f)
    if Logging.shouldlog(global_logger(), Logging.Debug, @__MODULE__, nothing, nothing)
        @eval import Plots
        display(Base.invokelatest(f))
    end
end

function parse_filter_scales(
    nmax::Integer,
    R0::Real;
    b::Real = √2
)
    @assert nmax ≥ 0 "nmax must be non-negative"
    nvals = 0:nmax
    return R0 .* (b .^ nvals)
end

function parse_filter_scales(
    scales::Union{Tuple{Vararg{Real}}, AbstractVector{<:Real}},
    R0::Real;   # R0 and b are unused, but kept for uniform call signature
    b::Real = √2
)
    return collect(scales)
end

@inline function apply_gaussian_filter!(
    buf::Array{ComplexF64,3},
    k2::Vector{Float64},
    R::Real,
)
    N = size(buf, 1)
    pref = -(R^2) / 2

    @inbounds for kz in 1:N, ky in 1:N, kx in 1:N
        ksq = k2[kx] + k2[ky] + k2[kz]
        buf[kx, ky, kz] *= exp(pref * ksq)
    end

    return buf
end

@inline function hessian_comp!(
    buf::Array{ComplexF64,3},
    tmp::Array{ComplexF64,3},
    k::Vector{Float64},
    dir1::Symbol,
    dir2::Symbol,
    R::Real,
)
    N = size(buf, 1)
    pref = -R^2

    @inbounds for iz in 1:N, iy in 1:N, ix in 1:N
        k1 = dir1 === :x ? k[ix] :
             dir1 === :y ? k[iy] :
                           k[iz]

        k2 = dir2 === :x ? k[ix] :
             dir2 === :y ? k[iy] :
                           k[iz]

        tmp[ix, iy, iz] = buf[ix, iy, iz] * (k1 * k2 * pref)
    end

    FFTW.ifft!(tmp)
    return tmp
end

@inline function copy_real_to_complex!(
    buf::AbstractArray{<:Complex}, 
    src::AbstractArray{<:Real}
    )
    @boundscheck size(buf) == size(src)
    @inbounds @simd for i in eachindex(buf, src)
        buf[i] = src[i]
    end
    return buf
end

@inline function copy_complex_to_real!(
    dst::AbstractArray{<:Real}, 
    buf::AbstractArray{<:Complex}
    )
    @boundscheck size(dst) == size(buf)
    @inbounds @simd for i in eachindex(dst, buf)
        dst[i] = real(buf[i])
    end
    return dst
end

@inline function compute_eigenvalues_sym3(
    a, 
    b, 
    c, 
    d, 
    e, 
    f)
    """Compute eigenvalues of a symmetric 3x3 matrix given its components and return in ascending order."""
    # matrix:
    # [ a  d  f
    #   d  b  e
    #   f  e  c ]

    p1 = d*d + e*e + f*f

    if p1 == 0.0
        # diagonal matrix
        l1 = a
        l2 = b
        l3 = c
    else
        q = (a + b + c) / 3
        a11 = a - q
        b22 = b - q
        c33 = c - q

        p2 = a11*a11 + b22*b22 + c33*c33 + 2*p1
        p  = sqrt(p2 / 6)

        # normalized matrix B = (1/p)(A - qI)
        r = (a11*b22*c33 + 2*d*e*f
             - a11*e*e - b22*f*f - c33*d*d) / (2*p^3)

        # guard for numerical error
        r = clamp(r, -1.0, 1.0)

        φ = acos(r) / 3

        l1 = q + 2*p*cos(φ)
        l3 = q + 2*p*cos(φ + 2π/3)
        l2 = 3*q - l1 - l3
    end

    # sort ascending
    if l1 > l2; l1, l2 = l2, l1; end
    if l2 > l3; l2, l3 = l3, l2; end
    if l1 > l2; l1, l2 = l2, l1; end

    return l1, l2, l3
end
