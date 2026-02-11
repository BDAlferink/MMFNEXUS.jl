function accumulate_components(
    labels, 
    ρ::Union{Nothing,AbstractArray}=nothing
    )
    maxid = maximum(labels)
    nt = Threads.maxthreadid()

    volume_t = [zeros(Float64, maxid) for _ in 1:nt]
    mass_t   = ρ === nothing ? nothing :
               [zeros(Float64, maxid) for _ in 1:nt]

    @threads for i in eachindex(labels)
        id = labels[i]
        id == 0 && continue

        tid = threadid()
        volume_t[tid][id] += 1
        if mass_t !== nothing
            mass_t[tid][id] += ρ[i]
        end
    end

    # Reduce
    volume = zeros(Float64, maxid)
    for t in 1:nt
        volume .+= volume_t[t]
    end

    if mass_t === nothing
        return volume, nothing
    end

    mass = zeros(Float64, maxid)
    for t in 1:nt
        mass .+= mass_t[t]
    end

    return volume, mass
end


valid_by_mass(volume, mass, min_mass) =
    findall(i -> volume[i] > 0 && mass[i] ≥ min_mass, eachindex(volume))

valid_by_volume(volume, min_vol) =
    findall(i -> volume[i] ≥ min_vol, eachindex(volume))

function mask_from_labels(
    labels::AbstractArray{<:Integer}, 
    valid_ids::Vector{Int}
    )
    maxid = maximum(labels)

    valid = falses(maxid)
    valid[valid_ids] .= true

    out = falses(size(labels))

    @inbounds for i in eachindex(labels)
        id = labels[i]
        out[i] = (id != 0) && valid[id]
    end

    return out
end

function node_properties(labels, 
    ρ::AbstractArray{<:Real,3},  
    sim_box::SimBox
    )
    volume, mass = accumulate_components(labels, ρ)

    volume .*= sim_box.VoxelVolume
    mass   .*= sim_box.VoxelVolume * sim_box.ρ_mean

    density_contrast = mass ./ volume ./ sim_box.ρ_mean
    return volume, mass, density_contrast
end

fraction_above_threshold(values, Δ) =
    sum(values .> Δ) / length(values)


"""
    node_threshold_fraction(sim_box::SimBox, ρ::Array{<:Real,3}, signature::Array{<:Real,3}, S_th::Real, min_node_mass::Real, Δ::Real)
    
    Compute the fraction of identified nodes above overdensity Δ at signature threshold S_th.
    This computation is used to find the optimal signature threshold for node identification.
"""
function node_threshold_fraction(
    sim_box, 
    ρ, 
    signature, 
    S_th, 
    min_node_mass, 
    Δ
)
    mask = signature .> S_th
    labels = label_components(mask)
    _ , mass, δ = node_properties(labels, ρ, sim_box)

    valid = mass .≥ min_node_mass
    return fraction_above_threshold(δ[valid], Δ)
end

"""
    insignificant_node_remover(sim_box::SimBox, ρ::Array{<:Real,3}, signature::Array{<:Real,3}, S_th::Real, min_node_mass::Real)
    
    Remove node objects below minimum mass from the signature field at given signature threshold S_th.
    Returns a boolean array indicating significant nodes.
"""
function insignificant_node_remover(
    sim_box::SimBox, 
    ρ::AbstractArray{<:Real,3}, 
    signature::AbstractArray{<:Real,3}, 
    S_th::Real, 
    min_node_mass::Real
    )
    # Label connected components above threshold
    mask = signature .> S_th
    labels = label_components(mask)

    volume, mass, _ = node_properties(labels, ρ, sim_box)

    # Only keep objects that have volume > 0 and a mass >= min_mass
    valid_ids = valid_by_mass(volume, mass, min_node_mass)
    valid_nodes = mask_from_labels(labels, valid_ids)

    return valid_nodes
end

function optimal_node_signature_threshold(
    sim_box, ρ, signatures, min_node_mass, Δ;
    s_low, s_high, fraction_target
)
    S_low  = maximum(signatures) * s_low
    S_high = maximum(signatures) * s_high

    # check if the given signature bounds are bounding the correct region
    if (node_threshold_fraction(sim_box, ρ, signatures, S_low, min_node_mass, Δ) > fraction_target) && (node_threshold_fraction(sim_box, ρ, signatures, S_high, min_node_mass, Δ) < fraction_target)
        error("Signature bounds are not valid")
        return nothing
    end

    f(x) = node_threshold_fraction(sim_box, ρ, signatures, x, min_node_mass, Δ) - fraction_target

    return fzero(f, (S_low, S_high), Roots.Brent(), atol=1e-3)
end


"""
    node_field(sim_box::SimBox, ρ::Array{<:Real,3}, signatures::Array{<:Real,3}, min_node_mass::Real, Δ::Real; s_low=1e-5, s_high=9e-1, fraction_target=0.5)

    Generate clean node field based on significant nodes above optimal signature threshold.
"""
function node_field(
    sim_box::SimBox, 
    ρ::AbstractArray{<:Real,3}, 
    signatures::AbstractArray{<:Real,3}, 
    min_node_mass::Real, 
    Δ::Real; 
    s_low=1e-5, 
    s_high=9e-1, 
    fraction_target=0.5
    )
    S_th = optimal_node_signature_threshold(
        sim_box, ρ, signatures, min_node_mass, Δ;
        s_low=s_low, s_high=s_high, fraction_target=fraction_target
    )

    labels = insignificant_node_remover(sim_box, ρ, signatures, S_th, min_node_mass)
    @info "Optimal node signature threshold: $(round(S_th, digits=3))"
    @info "Number of identified nodes: $(maximum(label_components(labels)))"
    clean_node = (signatures .* labels) .> S_th

    return clean_node
end

"""
    insignificant_fila_wall_remover(sim_box::SimBox, signatures::Array{<:Real,3}, S_th::Real, min_vol::Real)

    Remove filament/wall objects below minimum volume from the signature field at given signature threshold S_th.
    Returns a boolean array indicating significant filaments/walls.
"""
function insignificant_fila_wall_remover(
    sim_box::SimBox, 
    signatures::AbstractArray{<:Real,3}, 
    S_th::Real, 
    min_vol::Real
    )
    # Label connected components above threshold
    mask = signatures .> S_th
    labels = label_components(mask)

    volume, _ = accumulate_components(labels)

    volume .= Float32.(volume) .* sim_box.VoxelVolume

    # Only keep objects that have volume > min_vol
    valid_ids = valid_by_volume(volume, min_vol)
    valid_objects = mask_from_labels(labels, valid_ids)

    return valid_objects 
end

"""
    ΔMsq(ρ::Array{<:Real,3}, signatures::Array{<:Real,3}, S_range::Array{<:Real,1})

    Compute dM^2 / dlog10(S) for a range of signature thresholds S_range.
"""

using Base.Threads

function ΔMsq(
    ρ::AbstractArray{<:Real,3}, 
    signatures::AbstractArray{<:Real,3}, 
    S_range::AbstractArray{<:Real,1}
    )
    Mass = zeros(Float64, length(S_range))

    @threads for i in eachindex(S_range)
        S_th = S_range[i]
        s = 0.0
        @inbounds for j in eachindex(ρ)
            if signatures[j] > S_th
                s += ρ[j]
            end
        end
        Mass[i] = s
    end

    logS = log10.(S_range)
    dMsq = abs.(diff(Mass .^ 2) ./ diff(logS))
    midpoints = 10 .^ (logS[1:end-1] .+ diff(logS) ./ 2)

    return midpoints, dMsq
end

"""
    fila_wall_field(sim_box::SimBox, ρ::Array{<:Real,3}, signatures::Array{<:Real,3}, clean_field, min_vol::Real; s_low=1e-5, s_high=1e-2, stepsize=0.1)

    Generate clean filament/wall field based on significant filaments/walls above optimal signature threshold.
"""
function fila_wall_field(
    sim_box::SimBox, 
    ρ::AbstractArray{<:Real,3}, 
    signatures::AbstractArray{<:Real,3}, 
    clean_field, 
    min_vol::Real; 
    s_low=1e-5, 
    s_high=1e-2, 
    stepsize=0.1
    )
    # set signatures within already identified nodes and/or filaments to zero
    signatures[clean_field] .= 0
    
    S_low = maximum(signatures)*s_low
    S_high = maximum(signatures)*s_high

    # compute optimal signature threshold based on maximum dM^2/dlog10(S)
    signature_range = 10 .^ (log10.(S_low):stepsize:log10.(S_high))
    S_midpoints, Mass_diff_sq = ΔMsq(ρ, signatures, signature_range)

    # find maximum dM^2/dlog10(S) and corresponding signature threshold
    _ , idx = findmax(Mass_diff_sq)
    S_th = S_midpoints[idx]  
    
    # remove insignificant filaments/walls (too low volume) from the field
    fila_wall_labels = insignificant_fila_wall_remover(sim_box, signatures, S_th, min_vol)
    
    # Generate clean filament/wall field based on significant objects above the signature threshold
    clean_fila_wall = (signatures.*fila_wall_labels) .> S_th;

    return clean_fila_wall, S_th, S_midpoints, Mass_diff_sq./maximum(Mass_diff_sq)
end