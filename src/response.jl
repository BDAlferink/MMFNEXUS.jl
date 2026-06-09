"""
    FourierFilter(sim_box::SimBox)

Precomputes the k-grid and allocates reusable work arrays.
Call the object like a function to filter a density field in-place.
"""
struct FourierFilter

    k::Vector{Float64}
    k2::Vector{Float64}
    fftbuf::Array{ComplexF64,3}
    tmpbuf::Array{ComplexF64,3}
end

function FourierFilter(
    sim_box::SimBox
    )
    N = sim_box.N

    k = fftfreq(N) .* N .* 2π / sim_box.L
    k2 = k .^ 2
    
    fftbuf = Array{ComplexF64}(undef, N, N, N)
    tmpbuf = similar(fftbuf)

    FourierFilter(k, k2, fftbuf, tmpbuf)
end

"""
    (filt::FourierFilter)(ρ::Array{<:Real,3}, sim_box::SimBox, R::Real; log::Bool=false)

Apply the Gaussian Fourier filter to a real 3D density field ρ
with smoothing scale R. Reuses internal memory and performs FFT in-place.

Returns the filtered density field as a new array.
"""
function (filt::FourierFilter)(
    ρ::Array{<:Real,3}, 
    sim_box::SimBox, 
    R::Real, 
    log::Bool=false
    )
    N = sim_box.N
    buf = filt.fftbuf

    @assert size(ρ) == (N,N,N)

    # If log, work on a temporary copy so ρ remains intact
    work = log ? similar(ρ) : ρ

    # Take log10 if using NEXUS+ style filtering
    if log
        map!((x -> log10(x)), work, ρ)
    end
    
    # Copy real input into complex buffer
    copy_real_to_complex!(buf, work)

    # forward FFT
    FFTW.fft!(buf)

    # Gaussian filtering in Fourier space
    apply_gaussian_filter!(buf, filt.k2, R)
    
    # inverse FFT
    FFTW.ifft!(buf)

    # Output array
    out = similar(ρ)

    copy_complex_to_real!(out, buf)

    # exponentiate for NEXUS+ and renormalize by mean
    if log
        map!((x -> 10.0^x), out, out)
        out .*= (mean(ρ) / mean(out))
    end

    return out
end


"""
    fft_field!(filt::FourierFilter, ρ::Array{<:Real,3})

Fourier transform a real 3D field ρ, returning the Fourier-space field.
"""
function fft_field!(
    filt::FourierFilter, 
    ρ::Array{<:Real,3}
    )
    # Fourier transform a real 3D field ρ, returning the Fourier-space field.
    buf = filt.fftbuf
    copy_real_to_complex!(buf, ρ)
    FFTW.fft!(buf)
    return buf   # Returns Fourier-space field
end

"""
Compute signatures from Hessian eigenvalues.
"""
function compute_signatures!(S_n, 
    S_f, 
    S_w, 
    I, 
    l1, 
    l2, 
    l3
    )
    # precompute ratios for signatures
    invl1 = 1.0 / l1
    r21 = abs(l2 * invl1)
    r31 = abs(l3 * invl1)

    # negativity requirements
    m1 = ifelse(l1 < 0, 1.0, 0.0)
    m2 = ifelse(l2 < 0, 1.0, 0.0)
    m3 = ifelse(l3 < 0, 1.0, 0.0)
    c21 = ifelse(r21 < 1, 1.0, 0.0)
    c31 = ifelse(r31 < 1, 1.0, 0.0)

    # Only calculate required signatures

    # Node signature
    if S_n !== nothing
        S_n[I] = abs(l3*l3*invl1) * (m1*m2*m3)
    end

    # Filament signature
    if S_f !== nothing
        S_f[I] = abs(l2*l2*invl1) * (1 - r31) * (m1*m2*c31)
    end

    # Wall signature
    if S_w !== nothing
        S_w[I] = (1 - r21)*(1 - r31)*abs(l1) * (m1*c21*c31)
    end

    return nothing
end

function allocate_signature_arrays(
    N::Integer, 
    mode::Symbol
    )
    # Helper function to allocate only the signatures needed for NEXUS(+)
    if mode == :node
        return (
            S_n = Array{Float32}(undef, N, N, N),
            S_f = nothing,
            S_w = nothing
        )
    elseif mode == :fila_wall
        return (
            S_n = nothing,
            S_f = Array{Float32}(undef, N, N, N),
            S_w = Array{Float32}(undef, N, N, N)
        )
    else
        error("Unknown mode: $mode")
    end
end

struct HessianBuffers
    Hxx::Array{Float64,3}
    Hyy::Array{Float64,3}
    Hzz::Array{Float64,3}
    Hxy::Array{Float64,3}
    Hxz::Array{Float64,3}
    Hyz::Array{Float64,3}
end

function HessianBuffers(
    N::Int
    )
    HessianBuffers(
        Array{Float64}(undef, N,N,N),
        Array{Float64}(undef, N,N,N),
        Array{Float64}(undef, N,N,N),
        Array{Float64}(undef, N,N,N),
        Array{Float64}(undef, N,N,N),
        Array{Float64}(undef, N,N,N),
    )
end

"""
    signatures_hessian!(filt::FourierFilter, sim_box::SimBox, ρ::Array{<:Real,3}, R::Real; mode::Symbol)

Compute Hessian eigenvalue signatures (node, filament, wall) at smoothing scale R.
"""
function signatures_hessian!(
    filt::FourierFilter, 
    sim_box::SimBox, 
    ρ::Array{<:Real,3}, 
    R::Real, 
    hbuf::HessianBuffers; 
    mode::Symbol
    )

    N = sim_box.N
    # kx, ky, kz = filt.kx, filt.ky, filt.kz
    k = filt.k
    buf = filt.fftbuf
    tmp = filt.tmpbuf

    if mode == :fila_wall
        ρ_filtered = filt(ρ, sim_box, R, true)
    else
        ρ_filtered = ρ
    end

    # FFT the density field
    fft_field!(filt, ρ_filtered)

    # In case of NEXUS (node detection), we use a normal Gaussian filter. This can be done simultaneously
    # with the Hessian computation for efficiency. This is not possible with the logfilter.
    if mode == :node
        apply_gaussian_filter!(buf, filt.k2, R)
    end

    # Allocate only the signatures needed for NEXUS(+)
    S = allocate_signature_arrays(N, mode)

    S_n = S.S_n
    S_f = S.S_f
    S_w = S.S_w


    # Compute Hessian components into tmp, reuse storage each time. We immediately feed 
    # eigenvalue calculation voxel-by-voxel.

    hessian_comp!(buf, tmp, k, :x, :x, R)
    copy_complex_to_real!(hbuf.Hxx, tmp)

    hessian_comp!(buf, tmp, k, :y, :y, R)
    copy_complex_to_real!(hbuf.Hyy, tmp)

    hessian_comp!(buf, tmp, k, :z, :z, R)
    copy_complex_to_real!(hbuf.Hzz, tmp)

    hessian_comp!(buf, tmp, k, :x, :y, R)
    copy_complex_to_real!(hbuf.Hxy, tmp)

    hessian_comp!(buf, tmp, k, :x, :z, R)
    copy_complex_to_real!(hbuf.Hxz, tmp)

    hessian_comp!(buf, tmp, k, :y, :z, R)
    copy_complex_to_real!(hbuf.Hyz, tmp)

    # Compute eigenvalues + signatures voxel by voxel
    @inbounds for I in eachindex(hbuf.Hxx)
        l1, l2, l3 = compute_eigenvalues_sym3(hbuf.Hxx[I], hbuf.Hyy[I], hbuf.Hzz[I], hbuf.Hxy[I], hbuf.Hyz[I], hbuf.Hxz[I])

        compute_signatures!(S_n, S_f, S_w, I, l1, l2, l3)

    end

    return S_n, S_f, S_w
end

struct GaussKernels{K1, K2, K3}
    kernels::Tuple{K1, K2, K3}
end

function GaussKernels(R::Real)
    σ  = R
    hw = ceil(Int, 4σ)
    xs = -hw:hw

    g0 = exp.(.-xs.^2 ./ (2σ^2))
    g0 ./= cbrt(sum(g0)^3)

    mk(v, dim) = KernelFactors.ReshapedOneD{3, dim-1}(centered(v))

    kernels = ntuple(d -> mk(g0, d), 3)

    GaussKernels(kernels)
end


# ============================================================
# NEXUS  — nodes: Hessian of G_R * ρ
# ============================================================

"""
    compute_hessian_nodes!(hbuf, ρ, R; pad)

NEXUS node signature: H_ij = R² ∂ᵢ∂ⱼ (G_R * ρ)

deriv_method options:
    :finite_diff     — Smooth first, then 3-point finite-difference derivatives [default]
"""
function compute_hessian_nodes!(hbuf::HessianBuffers,
                                 ρ::AbstractArray{<:Real,3},
                                 R::Real;
                                 pad = "reflect")
    _compute_hessian_nodes_finite_diff!(hbuf, ρ, R; pad)
    return hbuf
end

# Smooth first, then finite-difference derivatives
function _compute_hessian_nodes_finite_diff!(hbuf::HessianBuffers,
                                              ρ::AbstractArray{<:Real,3},
                                              R::Real;
                                              pad = "reflect")
    gk = GaussKernels(R)
    ρ_R = imfilter(ρ, gk.kernels, pad)
    
    mk(v, dim) = KernelFactors.ReshapedOneD{3, dim-1}(v)
    d1f = ntuple(dim -> mk(_FD_D1, dim), 3)
    d2f = ntuple(dim -> mk(_FD_D2, dim), 3)
    idf = ntuple(dim -> mk(_FD_ID, dim), 3)
    
    R2 = R^2
    
    imfilter!(hbuf.Hxx, ρ_R, (d2f[1], idf[2], idf[3]), pad); hbuf.Hxx .*= R2
    imfilter!(hbuf.Hyy, ρ_R, (idf[1], d2f[2], idf[3]), pad); hbuf.Hyy .*= R2
    imfilter!(hbuf.Hzz, ρ_R, (idf[1], idf[2], d2f[3]), pad); hbuf.Hzz .*= R2
    imfilter!(hbuf.Hxy, ρ_R, (d1f[1], d1f[2], idf[3]), pad); hbuf.Hxy .*= R2
    imfilter!(hbuf.Hxz, ρ_R, (d1f[1], idf[2], d1f[3]), pad); hbuf.Hxz .*= R2
    imfilter!(hbuf.Hyz, ρ_R, (idf[1], d1f[2], d1f[3]), pad); hbuf.Hyz .*= R2
    
    return hbuf
end

# ============================================================
# NEXUS+ — filaments & walls: Hessian of f_R = C·10^(G_R * log10 ρ)
# ============================================================

"""
    compute_hessian_logfield!(hbuf, ρ, R; tmp_derivs, pad)

NEXUS+ filament/wall signature. Computes:

    g_R  = G_R * log10(ρ)
    f_R  = C · 10^g_R,   C chosen so ⟨f_R⟩ = ⟨ρ⟩

    H_ij = R² ∂ᵢ∂ⱼ f_R

"""
function compute_hessian_logfield!(hbuf::HessianBuffers,
                                    ρ::AbstractArray{<:Real,3},
                                    R::Real;
                                    tmp_derivs::Union{Nothing,NTuple{3,AbstractArray{<:Real,3}}} = nothing,
                                    pad = "reflect")
    _compute_hessian_logfield_finite_diff!(hbuf, ρ, R; tmp_derivs, pad)
    return hbuf
end

# Finite-difference derivatives on already-smoothed log field
function _compute_hessian_logfield_finite_diff!(hbuf::HessianBuffers,
                                                 ρ::AbstractArray{<:Real,3},
                                                 R::Real;
                                                 tmp_derivs::Union{Nothing,NTuple{3,AbstractArray{<:Real,3}}} = nothing,
                                                 pad = "reflect")
    # ln10 = log(10.0)
    # dims = size(ρ)
    gk = GaussKernels(R)
    
    f_R = exp10.(imfilter(log10.(ρ), gk.kernels, pad))
    f_R .*= mean(ρ) / mean(f_R)
        
    mk(v, dim) = KernelFactors.ReshapedOneD{3, dim-1}(v)
    d1f = ntuple(dim -> mk(_FD_D1, dim), 3)
    d2f = ntuple(dim -> mk(_FD_D2, dim), 3)
    idf = ntuple(dim -> mk(_FD_ID, dim), 3)
    
    R2 = R^2
    
    imfilter!(hbuf.Hxx, f_R, (d2f[1], idf[2], idf[3]), pad); hbuf.Hxx .*= R2
    imfilter!(hbuf.Hyy, f_R, (idf[1], d2f[2], idf[3]), pad); hbuf.Hyy .*= R2
    imfilter!(hbuf.Hzz, f_R, (idf[1], idf[2], d2f[3]), pad); hbuf.Hzz .*= R2
    imfilter!(hbuf.Hxy, f_R, (d1f[1], d1f[2], idf[3]), pad); hbuf.Hxy .*= R2
    imfilter!(hbuf.Hxz, f_R, (d1f[1], idf[2], d1f[3]), pad); hbuf.Hxz .*= R2
    imfilter!(hbuf.Hyz, f_R, (idf[1], d1f[2], d1f[3]), pad); hbuf.Hyz .*= R2

    return hbuf
end

"""
    signatures_hessian!(sim_box::SimBox, ρ::Array{<:Real,3}, R::Real; mode::Symbol)

Compute Hessian eigenvalue signatures (node, filament, wall) at smoothing scale R. Use a convolution in real space instead of Fourier filtering. This is provided for testing and comparison purposes.
"""
function signatures_hessian!(
    sim_box::SimBox, 
    ρ::Array{<:Real,3}, 
    R::Real, 
    hbuf::HessianBuffers; 
    mode::Symbol,
    pad = "reflect"
    )

    N = sim_box.N
    Δx = sim_box.L / sim_box.N
    R_grid = R / Δx

    if mode == :node
		compute_hessian_nodes!(hbuf, ρ, R_grid; pad=pad)
		
    elseif mode == :fila_wall
        compute_hessian_logfield!(hbuf, ρ, R_grid; pad=pad)
    else
        error("Unknown mode: $mode. Use :node or :fila_wall.")
    end

    # Allocate only the signatures needed for NEXUS(+)
    S = allocate_signature_arrays(N, mode)

    S_n = S.S_n
    S_f = S.S_f
    S_w = S.S_w

    # Compute eigenvalues + signatures voxel by voxel
    @inbounds for I in eachindex(hbuf.Hxx)
        l1, l2, l3 = compute_eigenvalues_sym3(hbuf.Hxx[I], hbuf.Hyy[I], hbuf.Hzz[I], hbuf.Hxy[I], hbuf.Hyz[I], hbuf.Hxz[I])

        compute_signatures!(S_n, S_f, S_w, I, l1, l2, l3)

    end

    return S_n, S_f, S_w
end

@inline function _update_max!(dst::Array{Float32,3}, src::Array{Float32,3})
    @inbounds @simd for i in eachindex(dst)
        v = src[i]
        if v > dst[i]; dst[i] = v; end
    end
end

"""
    multiscale_signature_max!(filt::FourierFilter, sim_box::SimBox, ρ::Array{<:Real,3}, filter_scales; mode::Symbol)

Compute maximum signatures across multiple scales.
"""
function multiscale_signature_max!(
    filt::FourierFilter, 
    sim_box::SimBox, 
    ρ::Array{<:Real,3}, 
    filter_scales; 
    mode::Symbol
    )
    N = sim_box.N

    S = allocate_signature_arrays(N, mode)

    # allocate output signature arrays
    S_n_max = S.S_n
    S_f_max = S.S_f
    S_w_max = S.S_w

    # Initialize max arrays to -Inf
    if S_n_max !== nothing
        fill!(S_n_max, -Inf32)
    end
    if S_f_max !== nothing
        fill!(S_f_max, -Inf32)
        fill!(S_w_max, -Inf32)
    end

    #allocate hessian buffers
    hbuf = HessianBuffers(N)

    # Loop over scales
    for R in filter_scales
        @info "Processing scale R = $(round(R, digits=3))"

        # compute signatures at scale R
        S_n, S_f, S_w = signatures_hessian!(filt, sim_box, ρ, R, hbuf; mode=mode)

        if mode == :node
            _update_max!(S_n_max, S_n)
        elseif mode == :fila_wall
            _update_max!(S_f_max, S_f)
            _update_max!(S_w_max, S_w)
        end
    end

    return S_n_max, S_f_max, S_w_max
end

"""
    multiscale_signature_max!(sim_box::SimBox, ρ::Array{<:Real,3}, filter_scales; mode::Symbol, pad)

Compute maximum signatures across multiple scales using real-space convolution.
Provided for comparison and testing; prefer the `FourierFilter` method for production use.
"""
function multiscale_signature_max!(
    sim_box::SimBox,
    ρ::Array{<:Real,3},
    filter_scales;
    mode::Symbol,
    pad = "reflect"
    )
    N = sim_box.N

    S = allocate_signature_arrays(N, mode)

    S_n_max = S.S_n
    S_f_max = S.S_f
    S_w_max = S.S_w

    if S_n_max !== nothing
        fill!(S_n_max, -Inf32)
    end
    if S_f_max !== nothing
        fill!(S_f_max, -Inf32)
        fill!(S_w_max, -Inf32)
    end

    hbuf = HessianBuffers(N)

    for R in filter_scales
        @info "Processing scale R = $(round(R, digits=3))"
        S_n, S_f, S_w = signatures_hessian!(sim_box, ρ, R, hbuf; mode=mode, pad=pad)

        if mode == :node
            _update_max!(S_n_max, S_n)
        elseif mode == :fila_wall
            _update_max!(S_f_max, S_f)
            _update_max!(S_w_max, S_w)
        end
    end

    return S_n_max, S_f_max, S_w_max
end