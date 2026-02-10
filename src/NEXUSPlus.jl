"""
    SimBox(N,L,M)

    Holds the simulation box parameters:
    N : Number of grid points per dimension
    L : Box size in cMpc/h
    M : Mass in Msun/h contained in the box
    Generates also additional useful parameters:
    ρ_mean : Mean density in Msun/h/(cMpc/h)^3
    VoxelVolume : Volume of a voxel in (cMpc/h)^3

"""
struct SimBox
    N::Int
    L::Float64
    M::Float64
    ρ_mean::Float64
    VoxelVolume::Float64

    function SimBox(
        N::Int, 
        L::Float64, 
        M::Float64
        )
        ρ_mean = M / L^3
        VoxelVolume = (L / N)^3
        new(N, L, M, ρ_mean, VoxelVolume)
    end
end


function NEXUS_Plus(
    ρ::AbstractArray{<:Real,3}, 
    N::Integer, 
    L::Real, 
    M::Real, 
    filter_parse, 
    Δ::Real, 
    min_node_mass::Real, 
    min_fila_volume::Real, 
    min_wall_volume::Real; 
    R0::Real = 0.5,
    level::Symbol = :info
)
    # Parse logging level. Options: 
    # - :debug elaborate with some figures
    # - :info standard output
    # - :none no output
    setup_logging(level)

    # Parse filter scales. Options:
    # - Integer: number of scales starting from R0 in sqrt(2) increments. R0 is default 0.5 Mpc/h but can be given as kwarg.
    # - Tuple/Array: user defined scales in Mpc/h
    filter_scales = parse_filter_scales(filter_parse, R0)

    # Initialize simulation box and Fourier filter
    sim_box = SimBox(N, L, M)
    filt = FourierFilter(sim_box)

    @info "Starting node signature calculation..."

    # find maximum signatures in scalespace for node identification, node option means NEXUS
    S_node_max, _, _ = multiscale_signature_max!(filt, sim_box, ρ, filter_scales; mode=:node);

    @info ""
    @info "Starting node identification..."
    # find node bool field given maximum signatures, contrast and min_mass
    # s_low and s_high are limiting threshold within which we look for the root. In case a root can not be found in this range the values could be changed. 
    # fraction_target is the optimal threshold for which half the objects are signficant as given in the original NEXUS paper
    clean_node = node_field(sim_box, ρ, S_node_max, min_node_mass, Δ; s_low=1e-5, s_high=9e-1, fraction_target=0.5);

    @info ""
    @info "Starting filament and wall signature calculation..."
    # find maximum signatures in scalespace for filament and wall identification, fila_wall option means NEXUS+ (logspace filter)
    _, S_fila_max, S_w_max = multiscale_signature_max!(filt, sim_box, ρ, filter_scales; mode=:fila_wall);

    @info ""
    @info "Starting filament identification..."
    # find filament bool field given maximum signatures and min_volume
    clean_fila, S_f_th, S_f_midpoints, f_Mass_diff_sq = fila_wall_field(sim_box, ρ, S_fila_max, clean_node, min_fila_volume; s_low=1e-5, s_high=9e-1, stepsize=0.1);
    @info "Optimal filament signature threshold: $(round(S_f_th, digits=3))"

    @info ""
    @info "Starting wall identification..."
    # find wall bool field given maximum signatures and min_volume
    clean_wall, S_w_th, S_w_midpoints, w_Mass_diff_sq = fila_wall_field(sim_box, ρ, S_w_max, clean_node .| clean_fila, min_wall_volume; s_low=1e-5, s_high=9e-1, stepsize=0.1);
    @info "Optimal wall signature threshold: $(round(S_w_th, digits=3))"


    @debug "Plotting ΔM² vs Signature"
    debugplot() do
        p = Plots.plot(
            S_f_midpoints,
            f_Mass_diff_sq;
            xscale = :log10,
            xlabel = "Signature",
            ylabel = "ΔM² (arbitrary units)",
            label = "Filaments",
            color = :blue
        )

        Plots.plot!(
            p,
            S_w_midpoints,
            w_Mass_diff_sq;
            label = "Walls",
            color = :green
        )

        display(p)
    end


    # find void bool field as the negation of all other environments
    clean_void = .!(clean_node .| clean_fila .| clean_wall);
    @info ""
    @info "NEXUS+ identification finished"

    @debug "Plot slice of the identified environments on top of density field"
    debugplot() do
        slice_index = div(sim_box.N, 2)

        x = range(0, sim_box.L; length = sim_box.N)
        y = range(0, sim_box.L; length = sim_box.N)

        p = Plots.heatmap(
            x,
            y,
            log10.(ρ[:, :, slice_index]);
            xlabel = "X [Mpc/h]",
            ylabel = "Y [Mpc/h]",
            color = :magma,
            clims = (minimum(log10.(ρ)), maximum(log10.(ρ))),
            xlims = (0, sim_box.L),
            ylims = (0, sim_box.L),
            aspect_ratio = :equal,
            label = "",  # No label for heatmap
            legend = :topright
        )
        Plots.contour!(
            p, x, y,
            clean_wall[:, :, slice_index];
            levels = [0.5],
            linewidth = 2,
            linecolor = :green,
            label = "Walls",
            colorbar_entry = false,
        )
        Plots.contour!(
            p, x, y,
            clean_fila[:, :, slice_index];
            levels = [0.5],
            linewidth = 2,
            linecolor = :blue,
            label = "Filaments",
            colorbar_entry = false,
        )
        Plots.contour!(
            p, x, y,
            clean_node[:, :, slice_index];
            levels = [0.5],
            linewidth = 2,
            linecolor = :red,
            label = "Nodes",
            colorbar_entry = false,
        )

            # Add invisible points for legend
        Plots.plot!(p, [], [], linewidth=2, linecolor=:green, label="Walls")
        Plots.plot!(p, [], [], linewidth=2, linecolor=:blue, label="Filaments")
        Plots.plot!(p, [], [], linewidth=2, linecolor=:red, label="Nodes")
        display(p)
    end


    return clean_node, clean_fila, clean_wall, clean_void
end