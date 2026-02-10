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

    return S_n, S_f, S_w
end

function allocate_signature_arrays(
    N::Integer, 
    mode::Symbol
    )
    # Helper functino to allocate only the signatures needed for NEXUS(+)
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

    # In case of NEXUS+, we use a log filter. This is used for filament and wall detection.
    # This does the filter of the logfield in fourier space before computing the hessian.
    if mode == :fila_wall
        gf = FourierFilter(sim_box)
        ρ_filtered = gf(ρ, sim_box, R, true)
    else
        ρ_filtered = ρ
    end

    # FFT the density field
    fft_field!(filt, ρ_filtered)

    # In case of NEXUS (node detection), we use a normal Gaussian filter. This can be done simultaniously
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

        S_n, S_f, S_w = compute_signatures!(S_n, S_f, S_w, I, l1, l2, l3)

    end

    return S_n, S_f, S_w
end

"""
    multiscale_signature_max!(filt::FourierFilter, sim_box::SimBox, ρ::Array{<:Real,3}, R0::Real; mode::Symbol, scalespec)
    Compute maximum signatures across multiple scales specified by the user.
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

        # update maximum signatures to the corresponding mode
        if mode == :node
            @inbounds @simd for i in eachindex(S_n_max)
                val = S_n[i]
                if val > S_n_max[i]
                    S_n_max[i] = val
                end
            end
        elseif mode == :fila_wall
            @inbounds @simd for i in eachindex(S_f_max)
                valf = S_f[i]
                if valf > S_f_max[i]
                    S_f_max[i] = valf
                end

                valw = S_w[i]
                if valw > S_w_max[i]
                    S_w_max[i] = valw
                end
            end
        end
    end

    return S_n_max, S_f_max, S_w_max
end