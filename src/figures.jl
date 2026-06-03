export genfig_tiny_io,
       genfig_results_rlf,
       genfig_results_intro_example_responses,
       genfig_results_intro_io_stack,
       genfig_results_intro_rlf_shift,
       genfig_results_cas_physio_Warren1989a_Fig6,
       genfig_results_cas_sim_tc,
       genfig_results_cas_sim_slf,
       genfig_results_cas_sim_slf_vs_params,
       genfig_results_cas_sim_innervation,
       genfig_results_cas_threshold_vs_cf,
       genfig_results_cas_threshold_vs_cf_supplemental,
       genfig_results_cas_ΔL_vs_cf,
       genfig_results_cas_ΔL_vs_cf_supplemental,
       genfig_results_cas_magnitude_vs_cf,
       genfig_summary_fit_to_physio,
       genfig_results_cas_sim_freq_sweeps,
       genfig_results_cas_sim_cf_vs_bsf,
       genfig_results_cas_sim_cf_vs_bsf_analyze,
       genfig_results_cas_band_widening,
       genfig_results_summary_rlf,
       genfig_maximum_attenuation,
       genfig_followup_elicitor_eps

# ##########################################################################################
# DEMO FIGURES
# Figures demonstrating aspects of the model, showing example model responses, etc.
# - genfig_tiny_io()
# - genfig_results_intro_example_responses()
# - genfig_results_intro_rlf()
# - genfig_results_intro_rlf_shift()
# - genfig_results_cas_sim_tc()
# - genfig_results_cas_sim_slf()
# - genfig_results_cas_sim_innervation()

"""
    genfig_tiny_io()

Plot tiny IO diagram for model
"""
function genfig_tiny_io(; fig=Figure(size=(180, 130)), ax=Axis(fig[1, 1]))
    # Plot MOC input-output nonlinearity example
    r = LinRange(0.0, 100.0, 500)
    γ = moc_nonlinearity.(r; w=1.0, β=0.05, minval=0.0, maxval=1.0, θ=1.0)
    lines!(ax, r, γ; color=:black, linewidth=4.0)
    hidedecorations!(ax)
    hidespines!(ax)

    # Display and save
    display(fig)
    saveplot("fig_tiny_io.png", fig)
    fig
end

"""
    genfig_results_intro_example_responses(params=MOCParams [; kwargs...])

Plot an example CAS stimulus waveform aset and responses w/ and w/o an elicitor
"""
function genfig_results_intro_example_responses(
    params::MOCParams=MOCParams3();
    lw=1.0, 
    config=Config(),
    fig=Figure(; size=(2*100, 4*100)),  # 2 in x 2 in
    level_elicitor=70.0,
    cf=8e3,
)
    # Stimulus (one without elicitor, one with elicitor)
    mp = midpoint(RLF(ipsi(fetch_model(params, cf; stage="hsr"))); config=config)  # estimate midpoint of ipsi RLF
    stims = map([-Inf, level_elicitor]) do lvl
         ProbePureToneElicitor(
            [cf, mp],  # place probe at midpoint
            [cf, lvl];  # place elicitor at -Inf or level_elicitor
            dur_probe=0.3, 
            dur_wait=0.05, 
            dur_delay=0.05, 
            dur_suppressor=0.2, 
            dur_post=0.05,
        )
    end

    # Plot stimuli (ipsi and contra)
    ax = Axis(fig[1, 1])
    xlims!(ax, 0.0, 0.40)
    ylims!(ax, -5e-4, 5e-4)
    hideydecorations!(ax)
    hidexdecorations!(ax, ticks=false)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(stims[1]), synthesize(ipsi(stims[2])); color=:black, linewidth=lw)

    ax = Axis(fig[2, 1])
    xlims!(ax, 0.0, 0.40)
    ylims!(ax, -15e-2, 15e-2)
    hideydecorations!(ax)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(stims[1]), synthesize(contra(stims[2])); color=:black, linewidth=lw)

    # Create LPF to smooth lines
    lpf = digitalfilter(Lowpass(0.5e3; fs=100e3), Butterworth(4))

    # Plot response w/o elicitor
    ax = Axis(fig[3, 1])
    r = @memo config compute(fetch_model(params, cf; stage="hsr"), stims[1])
    ylims!(ax, 0.0, 1000.0)
    xlims!(ax, 0.0, 0.40)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(r, samprate(stims[1])), filtfilt(lpf, r), color=:gray, linewidth=lw/2, linestyle=:dash)

    # Plot response w/ elicitor
    r = @memo config compute(fetch_model(params, cf; stage="hsr"), stims[2])
    ylims!(ax, 0.0, 400.0)
    xlims!(ax, 0.0, 0.40)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(r, samprate(stims[2])), filtfilt(lpf, r), color=:red, linewidth=lw)

    # Plot MOC response
    ax = Axis(fig[4, 1]; xticklabelsvisible=false)
    r = @memo config Utilities._compute(contra(fetch_model(params, cf)), contra(stims[2]))
    r = r["mocwdr"][contra(fetch_model(params, cf)).coi[1]]  # extract MOC WDR gain control signal
    ylims!(ax, 0.0, 80.0)
    xlims!(ax, 0.0, 0.40)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(r, samprate(stims[2])), r, color=:black, linewidth=lw)

    # Plot gain factor
    ax = Axis(fig[5, 1]; xticklabelsvisible=false)
    r = @memo config compute(fetch_model(params, cf; stage="gainpostmix"), stims[2])
    ylims!(ax, 0.0, 1.25)
    xlims!(ax, 0.0, 0.40)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    ax.yticks = [0.0, 0.5, 1.0]
    lines!(ax, timevec(r, samprate(stims[2])), r, color=:black, linewidth=lw)

    # Plot gain dB
    ax = Axis(fig[6, 1])
    ylims!(ax, 0.0, 10.0)
    xlims!(ax, 0.0, 0.40)
    vlines!(ax, [0.1]; color=:gray, linewidth=lw/2)
    lines!(ax, timevec(r, samprate(stims[2])), cohc_to_ΔL(r), color=:black, linewidth=lw)

    # Adjust row/gap size
    rowsize!(fig.layout, 1, Relative(0.15))  # stim row
    rowsize!(fig.layout, 2, Relative(0.15))  # stim row
    rowgap!(fig.layout, Relative(0.02))
    rowgap!(fig.layout, 2, Relative(0.05))
    rowgap!(fig.layout, 3, Relative(0.05))

    # Display and return
    display(fig)
    saveplot("fig_results_intro_response_examples.png", fig)
    fig
end

function genfig_results_intro_io_stack(
    params::MOCParams,
    params_unweighted::MOCParams4;
    lw=1.0, 
    config=Config(),
    fig=Figure(; size=(1.8*100, 4*100)),  # 2 in x 2 in
    level_elicitor=65.0,
    cf=8e3,
)
    # First begin by plotting rate-to-γ nonlinearity
    ax = Axis(fig[1, 1], xlabel="MOC rate (sp/s)", ylabel="Gain factor (γ)")
    r = LinRange(0.0, 100.0, 500)
    γ = moc_nonlinearity.(
        r; 
        w=params.moc_weight(cf), 
        β=params.moc_beta(cf), 
        minval=0.0, 
        maxval=1.0, 
        θ=params.moc_offset(cf)*params.moc_weight(cf),
    )
    map([40.0, 60.0, 80.0]) do lvl
        stim = PureTone(; freq=cf, level=lvl, dur=1.0)
        model = fetch_model_baseline(cf; stage="mocwdr")
        resp = compute(model, stim)
        m = maximum(resp)
        vlines!(ax, [m]; color=:gray, linewidth=1.0)
    end
    lines!(ax, r, γ; color=:black)

    # Next, plot empirical ΔL vs γ
    ax = Axis(fig[2, 1], xlabel="Gain factor (γ)", ylabel="ΔL (dB)")
    γ = LinRange(0.0, 1.0, 500)
    ΔL = cohc_to_ΔL(γ; cf=cf, config=config)
    lines!(ax, γ, ΔL; color=:black)

    # Plot variation between offset and CF
    ax = tradlogxax(fig[3, 1]; ylabel="Offset / θ (sp/s)", xticklabelsvisible=false, xlabelvisible=false)
    cfs = LogRange(0.25e3, 20e3, 500)
    lines!(ax, cfs ./ 1e3, params.moc_offset.(cfs); color=:black)
    lines!(ax, cfs ./ 1e3, params_unweighted.moc_offset.(cfs); color=:gray, linestyle=:dash)

    # Plot variation between β and CF
    ax = tradlogxax(fig[4, 1]; xlabel="CF (kHz)", ylabel="Slope / β")
    cfs = LogRange(0.25e3, 20e3, 500)
    lines!(ax, cfs ./ 1e3, params.moc_beta.(cfs); color=:black)
    lines!(ax, cfs ./ 1e3, params_unweighted.moc_beta.(cfs); color=:gray, linestyle=:dash)

    # Adjust gaps
    rowgap!(fig.layout, 1, Relative(0.04))
    rowgap!(fig.layout, 2, Relative(0.04))
    rowgap!(fig.layout, 3, Relative(0.02))

    # Display and return
    display(fig)
    saveplot("fig_results_intro_io_stack.png", fig)
    fig
end

"""
    genfig_results_rlf()

Plot example rate-level function with annotations for threshold and midpoint
"""
function genfig_results_rlf(; config=Config())
    # Set up simulation
    model = fetch_model_baseline(2e3)
    sim = RLF(model; levels=-10.0:2.0:80.0)

    # Create figure
    fig = Figure(; size=(2*100, 1.7*100))
    ax = Axis(fig[1, 1])

    # Simulate RLF
    μ = @memo config simulate(sim)

    # Smooth rlf
    l̂, μ̂ = smooth(sim, level(sim), μ)

    # Add vertical lines at threshold and midpoint
    θ = Utilities.threshold(sim; config=config, interp=true)
    mp = Utilities.midpoint(sim; config=config, interp=true)
    idx_θ = argmin(abs.(l̂ .- θ))
    idx_mp = argmin(abs.(l̂ .- mp))
    lines!(ax, [-20.0, θ], [μ̂[idx_θ], μ̂[idx_θ]]; color=:gray, linewidth=1.0)
    arrows!(ax, [θ], [μ̂[idx_θ]], [0.0], [-μ̂[idx_θ]+20.0]; color=:gray, linewidth=1.0)
    lines!(ax, [-20.0, mp], [μ̂[idx_mp], μ̂[idx_mp]]; color=:gray, linewidth=1.0)
    arrows!(ax, [mp], [μ̂[idx_mp]], [0.0], [-μ̂[idx_mp]+20.0]; color=:gray, linewidth=1.0)

    # Plot RLF over arrows
    lines!(ax, l̂, μ̂; color=:black)

    # Calculate DR and plot
    dr, levels, idxs = dynamicrange(sim; config=config, interp=true)
    rates = (μ̂[idxs[1]], μ̂[idxs[2]]) 
    x = 35.0
    lines!(ax, [x, x+2], [rates[1], rates[1]]; color=:gray)
    lines!(ax, [x, x+2], [rates[2], rates[2]]; color=:gray)
    lines!(ax, [x+2, x+2], [rates[1], rates[2]]; color=:gray)

    y = 340.0
    lines!(ax, [levels[1], levels[1]], [y, y-10.0]; color=:gray)
    lines!(ax, [levels[2], levels[2]], [y, y-10.0]; color=:gray)
    lines!(ax, [levels[1], levels[2]], [y, y]; color=:gray)

    # Scale
    ylims!(ax, 0.0, 360.0)
    xlims!(ax, -5.0, 45.0)

    # Add labels
    ax.xlabel = "Probe level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"

    # Plot and save
    display(fig)
    saveplot("fig_results_intro_rlf_clean.png", fig)
    fig
end

"""
    genfig_results_intro_rlf_shift()

Plot example of a single RLF shift with COHC=0.4
"""
function genfig_results_intro_rlf_shift(
    params::MOCParams=MOCParams3(); 
    config=Config(),
    level_elicitor=65.0,
)
    # Set up simulation
    model = fetch_model(params, 2e3)
    sim_noeff = RLF(model, ProbePureToneElicitor3; levels=-10.0:2.0:80.0, level_elicitor=-Inf)
    sim_witheff = RLF(model, ProbePureToneElicitor3; levels=-10.0:2.0:80.0, level_elicitor=level_elicitor)

    # Create figure
    fig = Figure(; size=(2*100, 1.7*100))  # 2 in x 1.7 in
    ax = Axis(fig[1, 1])

    # Extract thresholds and midpoints
    th_noeff, th_witheff = Utilities.threshold.([sim_noeff, sim_witheff]; config=config)
    mp_noeff, mp_witheff = midpoint.([sim_noeff, sim_witheff]; config=config)

    # Also grab rates while we're at it
    μ_noeff = @memo config simulate(sim_noeff)
    μ_witheff = @memo config simulate(sim_witheff)
    l = level(sim_noeff)

    # Plot shifts
    r_mp_noeff = μ_noeff[argmin(abs.(mp_noeff .- level(sim_noeff)))]
    r_mp_witheff = μ_witheff[argmin(abs.(mp_noeff .- level(sim_witheff)))]
    arrows!(ax, [mp_noeff], [r_mp_noeff], [mp_witheff - mp_noeff - 2.0], [0.0]; color=:gray)
    arrows!(ax, [mp_noeff], [r_mp_noeff], [0.0], [r_mp_witheff - r_mp_noeff + 20.0]; color=:gray)

    # Plot RLFs
    lines!(ax, level(sim_noeff), μ_noeff; color=:black)
    lines!(ax, level(sim_witheff), μ_witheff; color=:red)

    # Scale
    ylims!(ax, 0.0, 360.0)
    xlims!(ax, -5.0, 45.0)

    # Add labels
    ax.xlabel = "Probe level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"

    # Plot and save
    display(fig)
    saveplot("fig_results_intro_rlf_shift.png", fig)
    fig
end

"""
    genfig_results_cas_sim_tc()

Plot example tuning curve from afferent HSR fiber.

This function is used to generate Figure 3A and Figure 4A. Here, an iso-response tuning
curve is estimated by collecting rate-level functions at several frequencies around the CF
and then estimating a threshold from each using the established threshold-estimation
strategies used for RLFs. The CF is fixed at 2 kHz and the model type is fixed at HSR.
"""
function genfig_results_cas_sim_tc(; 
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogxax(fig[1, 1]),
    fn="fig_results_cas_sim_tuning_curve.png",
    annotate_flag=1,
    level_high=100.0,
    stage="hsr",
    xlims=(0.1, 10.0),
    savefig=false,
)
    # Fetch 2 kHz model
    model = fetch_model(2e3; stage=stage)

    # Set up RLF tuning curve simulation and extract threshold curve from it, plot into axis
    sim = RLFTC(ipsi(model); level_high=level_high, freq_low=-3.0, freq_high=3.0, n_freq=71)
    f = Utilities.freq(sim)
    θ = threshold_curve(sim; config=config)
    x̂, ŷ = smooth(log2.(f), θ, :loess, span=0.15, degree=2)
    lines!(ax, (2 .^ x̂) ./ 1e3, ŷ; color=:black, linewidth=1.0)

    # Optionally, annotate different stimulus paradigms related to CAS
    @match annotate_flag begin
        1 => begin
            # Case 1: Annotate SLF paradigm (mark fixed probe and level-varying elicitor)
            mp = midpoint(RLF(ipsi(model); levels=0.0:1.0:55.0); config=config)
            scatter!(ax, [getcf(model)/1e3], [mp]; color=:gray, markersize=5.0)
            lines!(ax, [getcf(model)/1e3, getcf(model)/1e3], [20.0, 60.0]; color=:red)
        end
        2 => begin
            # Case 2: Annotate ??
            lines!(ax, [getcf(model)/1e3, getcf(model)/1e3], [0.0, 80.0]; color=:gray)
            scatter!(ax, [getcf(model)/1e3], [60.0]; color=:red, markersize=5.0)
        end
        3 => begin
            # Case 3: Annotate frequency sweep paradigm (mark fixed probe and freq-varying elicitor)
            mp = midpoint(RLF(ipsi(model); levels=0.0:1.0:55.0); config=config)
            scatter!(ax, [getcf(model)/1e3], [mp]; color=:gray, markersize=5.0)
            scatter!(ax, [getcf(model)/1e3 / 4.0, getcf(model)/1e3 * 4.0], [60.0, 60.0]; color=:red)
            lines!(ax, [getcf(model)/1e3 / 4.0, getcf(model)/1e3 * 4.0], [60.0, 60.0]; color=:red)
        end
    end

    # Adjust scale and ticks
    xlims!(ax, xlims...)
    # ax.xminorticks = [0.5, 0.6, 0.7, 0.8, 0.9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    ax.xticks = [0.1, 1.0, 10.0]
    ylims!(ax, -10.0, 90.0)

    # Add labels
    ax.xlabel = "Frequency (kHz)"
    ax.ylabel = "Threshold (dB SPL)"

    # Display and save
    display(fig)
    if savefig
        saveplot(fn, fig)
    end
    fig
end

"""
    genfig_results_cas_sim_slf()

Plot example of suppressor-level function.

Plots example suppressor-level function for 2-kHz ipsi CF in the HSR model. The ipsi probe
is presented at the midpoint of the ipsi RLF (with 2-dB resolution) on CF. The contra probe
is presented on CF with a wide level range.
"""
function genfig_results_cas_sim_slf(
    params::MOCParams=MOCParams3();
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    savefig=false,
)
    # Plot
    model = fetch_model(params, 2e3)
    sim = SLF(
        model; 
        level_probe=midpoint(RLF(ipsi(model); levels=0.0:2.0:60.0); config=config), 
        level_elicitor_min=15.0,
        level_elicitor_max=80.0,
    )
    r = @memo config simulate(sim)
    scatter!(ax, level(sim), r; color=:black, markersize=4.0)
    lines!(ax, level(sim), r; color=:black)

    # Scale
    ylims!(ax, 0.0, 250.0)
    xlims!(ax, 10.0, 85.0)

    # Add labels
    ax.xlabel = "Elicitor level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"

    # Add indicator for threshold
    th = threshold(sim; config=config)
    arrows!(ax, [th], [250.0], [0.0], [-25.0]; color=:red)

    # Display and save
    display(fig)
    if savefig
        saveplot("fig_results_cas_sim_slf.png", fig)
    end
end

function genfig_results_cas_sim_slf_vs_params(
    params::MOCParams=MOCParams3();
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    β_mults=[1.0, 1.5, 3.0],
    θ_adds=[0.0, 2.0, 4.0],
    colors=colorschemes[:Dark2_8],
    savefig=false,
)
    # Plot
    map(zip(β_mults, colors)) do (β_mult, color)
        params_adjusted = MOCParams4(;
            moc_width=params.moc_width,
            moc_beta=x -> params.moc_beta(x) * β_mult,
            moc_offset=params.moc_offset,
            moc_weight=params.moc_weight,
        )
        model = fetch_model(params_adjusted, 2e3)
        sim = SLF(
            model; 
            level_probe=midpoint(RLF(ipsi(model); levels=0.0:2.0:60.0); config=config), 
            level_elicitor_min=15.0,
            level_elicitor_max=80.0,
        )
        r = @memo config simulate(sim)
        scatter!(ax, level(sim), r; color=color, markersize=4.0)
        lines!(ax, level(sim), r; color=color)
        th = threshold(sim; config=config)
        r_at_th = r[argmin(abs.(th .- level(sim)))]
        arrows!(ax, [th], [250.0], [0.0], [-25.0]; color=color)
    end

    # Scale
    ylims!(ax, 0.0, 250.0)
    xlims!(ax, 10.0, 85.0)

    # Add labels
    ax.xlabel = "Elicitor level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"

    # Add indicator for threshold

    # Display and save
    display(fig)
    if savefig
        saveplot("fig_results_cas_sim_slf_vs_params.png", fig)
    end
end


"""
    genfig_results_cas_sim_innervation(params::MOCParams[; size=(2.0*100, 1.6*100), fig=Figure(), ax=Axis()])

Plot MOC weight vs CF 
"""
function genfig_results_cas_sim_innervation(
    params::Vector{<:MOCParams}=[MOCParams3()];
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogxax(fig[1, 1]),
    linestyles=[:dash, :solid],
)
    # Select CF axis and plot weight
    cf = LogRange(0.01e3, 32e3, 500)
    map(zip(params, linestyles)) do (p, l)
        lines!(ax, cf ./ 1e3, p.moc_weight_func.(cf); color=:black, linestyle=l)
    end

    # Adjust ticks
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "MOC weight"

    # Adjust limits
    ylims!(ax, 0.0, 1.1*maximum(params[1].moc_weight_func.(cf)))
    xlims!(ax, 0.1, 30.0)

    # Display and save
    display(fig)
    saveplot("fig_results_cas_sim_innervation.png", fig)
end

function genfig_results_cas_sim_innervation(params::MOCParams, args...; kwargs...)
    genfig_results_cas_sim_innervation([params], args...; kwargs...)
end


# ##########################################################################################
# PHYSIO FIGURES
# Figures plotting physio plots either de novo or via wrappers defined in `data.jl` files
# - genfig_results_cas_physio_Warren1989a_Fig6()

"""
    genfig_results_cas_physio_Warren1989a_Fig6()

Plot Figure 6 from Warren and Liberman (1989a) as reference data
"""
function genfig_results_cas_physio_Warren1989a_Fig6()
    # Create figure 
    fig = Figure(; size=(2.5*100, 2*100))  # 2.5" x 2"

    # Plot physiological data as base layer
    _, ax = plot_Warren1989a_Fig6(; 
        fig=fig, 
        markersize=8.0,
        linewidth=2.5,
    )
    xlims!(ax, 0.2, 40.0)
    ylims!(ax, 10.0, 90.0)

    # Display and save
    display(fig)
    saveplot("fig_results_cas_physio_Warren1989a_Fig6.png", fig)
end


# ##########################################################################################
# SIMULATION FIGURES
# Figures showing simulation results. These constitute the primary figures in the paper
# and they all use a similar syntax.
# - genfig_results_cas_threshold_vs_cf()
# - genfig_results_cas_ΔL_vs_cf()
# - genfig_results_cas_magnitude_vs_cf()
# - genfig_results_cas_sim_freq_sweeps()
# - genfig_results_cas_sim_cf_vs_bsf()


"""
    genfig_results_cas_threshold_vs_cf([params=MOCParams(); config=Config()])

Plot CAS threshold vs CF in comparison to physio data from Warren and Liberman (1989a)

Plot CAS thresholds for HSR model in comparison to simplified visualization of Warren and
Liberman's CAS threshold data. As a useful reference, we also plot the absolute threshold 
for the LSR model at the same CFs; this makes it easy to see how CAS thresholds closely
match afferent thresholds if no weighting function is applied.
"""
function genfig_results_cas_threshold_vs_cf(
    params::MOCParams; 
    # Configuration
    config=Config(),
    # Stimulus/simulation parameters
    n_cf=21,
    level_step_elicitor=2.5,
    interp=true,
    # Figure
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    color_physio=(fetch_attr(:color, :physio), 0.2),
    linewidth=2.5,
    color=:black,
    plot_midpoint=false,
    plot_physio=true,
    marker=:circle,
    markersize=4.5,
    savefig=false,
    ylims=(10.0, 90.0),
)
    # Plot physiological data as base layer
    if plot_physio
        plot_Warren1989a_Fig6(; 
            fig=fig,
            ax=ax,
            markersize=6.0,
            linewidth=2.5,
            colors=fill(color_physio, 3),
        )
    end

    # Plot absolute threshold across relevant CF range for LSR baseline model 
    cfs = LogRange(0.5e3, 16e3, n_cf)
    θ = lsr_threshold.(cfs; config=config)
    lines!(ax, cfs ./ 1e3, θ; color=fetch_attr(:color, :lsr), linestyle=:dot, linewidth=linewidth/2)

    # Optionally plot midpoint across relevant CF range for HSR baseline model (this is where we place the ipsi probe)
    if plot_midpoint
        cfs = LogRange(0.5e3, 16e3, n_cf)
        mp = hsr_midpoint.(cfs; config=config)
        lines!(ax, cfs ./ 1e3, mp; color=fetch_attr(:color, :hsr), linestyle=:dot, linewidth=linewidth/2)
    end

    # Pick CFs ranging from 0.5 to 16 kHz, spaced evenly on log range
    cfs = LogRange(0.5e3, 16e3, n_cf)

    # Select models
    models = map(x -> fetch_model(params, x; stage="hsr"), cfs)

    # Run SLF for each model, placing probe at midpoint of ipsilateral RLF
    sims = map(models) do model
        SLF(
            model; 
            levels_elicitor=10.0:level_step_elicitor:70.0, 
            level_probe=hsr_midpoint(getcf(model); config=config),
        )
    end

    # Extract thresholds and plot
    θ = Utilities.threshold.(sims; config=config, interp=interp)
    scatter!(ax, cfs ./ 1e3, θ; color=color, markersize=markersize, marker=marker)

    # Adjust xlims and ylims to leave more room for action
    ylims!(ax, ylims...)
    xlims!(ax, 0.2, 30.0)
    ax.yticks = 20.0:10.0:90.0

    # Display and save
    display(fig)
    if savefig
        saveplot("fig_results_cas_sim_thresholds.png", fig) 
    end
    return fig, ax, sims
end

function genfig_results_cas_threshold_vs_cf(
    params::Vector{<:MOCParams}, 
    args...; 
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    markers=[:cross, :circle],
    colors=[(fetch_attr(:color, :model), 0.5), fetch_attr(:color, :model)],
    kwargs...
)
    for (idx, (p, m, c)) in enumerate(zip(params, markers, colors))
        genfig_results_cas_threshold_vs_cf(p; marker=m, color=c, fig=fig, ax=ax, plot_physio=idx==1, kwargs...)
    end
end

function genfig_results_cas_threshold_vs_cf_supplemental(
    sims;
    config=Config(),
)
    # Create a figure with one small panel for each CF
    fig = Figure(; size=(150*length(sims), 200))

    # Loop over sims and plot
    for (idx, sim) in enumerate(sims)
        ax = Axis(fig[1, idx]; yticklabelsvisible=idx==1, title="CF = $(round(getcf(sim)/1e3; sigdigits=2)) kHz")
        viz(sim; config=config, fig=fig, ax=ax)
    end
    fig
end

"""
    genfig_results_cas_ΔL_vs_cf([params=MOCParams3(); config=Config()])

Plot CAS ΔL vs CF in comparison to physio data from Warren and Liberman (1989b)

Plot CAS ΔLs for HSR model in comparison to simplified visualization of Warren and
Liberman's CAS threshold data.
"""
function genfig_results_cas_ΔL_vs_cf(
    params::MOCParams;
    config=Config(),
    level_step_probe=2.5,
    n_cf=21,
    level_elicitor=65.0,
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]; xscale=log10),
    marker=:circle,
    color=:black,
    savefig=false,
    kwargs...
)
    # Plot physiological data as base layer
    plot_Warren1989b_Fig2(; fig=fig, ax=ax, markersize=4.0, linewidth=2.5, color=:gray)

    # Set limits
    xlims!(ax, 0.1, 30.0)
    ylims!(ax, 0.0, 20.0)

    # Pick CFs ranging from 0.5 to 16 kHz, spaced evenly on log range
    cfs = LogRange(0.5e3, 16e3, n_cf)

    # Select models
    models = map(x -> fetch_model(params, x; stage="hsr"), cfs)

    # For each CF/model, measure ΔL at fixed elicitor level relative to no elicitor
    sims = map(models) do model
        rlf_baseline = RLF(
            model,
            ProbePureToneElicitor3; 
            level_step=level_step_probe,
            level_elicitor=-Inf,
        )

        rlf_target = RLF(
            model,
            ProbePureToneElicitor3; 
            level_step=level_step_probe,
            level_elicitor=level_elicitor,
        )
        return rlf_baseline, rlf_target
    end
    
    # Get all ΔL for plotting
    ΔL = map(x -> midpoint_shift(x[1], x[2]; config=config, interp=true), sims)

    # Also get average dynamic ranges and rate ranges
    dr_without = map(x -> dynamicrange(x[1]; config=config, interp=true)[1], sims)
    dr_with = map(x -> dynamicrange(x[2]; config=config, interp=true)[1], sims)
    rr_without = map(x -> raterange(x[1]; config=config), sims)
    rr_with = map(x -> raterange(x[2]; config=config), sims)
    println("(DR, RR) without = ($(mean(skipnan(dr_without))), $(mean(rr_without)))")
    println("(DR, RR) with = ($(mean(skipnan(dr_with))), $(mean(rr_with)))")

    # Extract RLF shifts and plot
    scatter!(ax, cfs ./ 1e3, ΔL; color=color, markersize=4.5, marker=marker)

    # Display and save
    display(fig)
    if savefig 
        saveplot("fig_results_cas_sim_shifts.png", fig) 
    end
    fig, ax, sims
end

function genfig_results_cas_ΔL_vs_cf(
    params::Vector{<:MOCParams}, 
    args...; 
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    markers=[:cross, :circle],
    colors=[(fetch_attr(:color, :model), 0.5), fetch_attr(:color, :model)],
    kwargs...
)
    for (idx, (p, m, c)) in enumerate(zip(params, markers, colors))
        genfig_results_cas_ΔL_vs_cf(p; marker=m, color=c, fig=fig, ax=ax, plot_physio=idx==1, kwargs...)
    end
end

function genfig_results_cas_ΔL_vs_cf_supplemental(
    sims;
    config=Config(),
)
    # Create a figure with one small panel for each CF
    fig = Figure(; size=(150*length(sims), 200))

    # Loop over sims and plot
    for (idx, (rlf_baseline, rlf_precursor)) in enumerate(sims)
        ax = Axis(fig[1, idx]; yticklabelsvisible=idx==1, title="CF = $(round(getcf(rlf_baseline)/1e3; sigdigits=2)) kHz", ylabelvisible=idx==1)
        viz(rlf_baseline; config=config, fig=fig, ax=ax, color=:black, markersize=6.0, interp=true)
        viz(rlf_precursor; config=config, fig=fig, ax=ax, color=:red, markersize=6.0, interp=true)
    end
    fig
end

"""
    genfig_results_cas_magnitude_vs_cf([params=MOCParams3(); config=Config()])

Plot CAS magnitude vs CF in comparison to physio data from Warren and Liberman (1989a)

Plot CAS magnitudes for HSR model in comparison to simplified visualization of Warren and
Liberman's CAS magnitude data. Magnitude is calculated at a constant supra-threshold 
elicitor level of 65 dB SPL.
"""
function genfig_results_cas_magnitude_vs_cf(
    params::MOCParams;
    config=Config(),
    n_cf=21,
    level_elicitor=70.0,
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    plot_physio=true,
    savefig=false,
    color=:black,
    linestyle=:solid,
    marker=:circle,
    kwargs...
)
    # Plot physiological data as base layer
    if plot_physio
        plot_Warren1989a_Fig7(; 
            fig=fig,
            ax=ax,
            markersize=6.0,
            linewidth=2.5,
            colors=[:gray, :gray, :gray],
        )
    end

    # Plot suppression magnitude either WITHOUT or WITH CF weighting
    # Pick CFs ranging from 0.5 to 16 kHz, spaced evenly on log range
    cfs = LogRange(0.5e3, 16e3, n_cf)

    # Select models
    models = map(x -> fetch_model(params, x; stage="hsr"), cfs)

    # Run SLF for each model, placing probe at midpoint of ipsilateral RLF
    sims = map(models) do model
        SLF(
            model; 
            levels_elicitor=[-Inf, level_elicitor], 
            level_probe=hsr_midpoint(getcf(model); config=config)
        )
    end

    # Extract thresholds and plot
    m = magnitude.(sims; config=config)
    scatter!(ax, cfs ./ 1e3, m; color=color, markersize=4.5, marker=marker)
    # x̂, ŷ = smooth(log2.(cfs), m, quadratic)
    # lines!(ax, 2 .^ x̂ ./ 1e3, ŷ; color=color, linewidth=2.5, linestyle=linestyle)

    ylims!(ax, 100.0, -20.0)

    # Display and save
    display(fig)
    if savefig saveplot("fig_results_cas_sim_magnitudes.png", fig) end
    return fig, ax, sims
end

function genfig_results_cas_magnitude_vs_cf(
    params::Vector{<:MOCParams}, 
    args...; 
    linestyles=[:dash, :solid],
    markers=[:cross, :circle],
    colors=[(fetch_attr(:color, :model), 0.5), fetch_attr(:color, :model)],
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    kwargs...
)
    for (idx, (p, l, m, c)) in enumerate(zip(params, linestyles, markers, colors))
        genfig_results_cas_magnitude_vs_cf(p, args...; linestyle=l, marker=m, color=c, fig=fig, ax=ax, plot_physio=idx==1, kwargs...)
    end
end


"""
    genfig_summary_fit_to_physio([params=MOCParams3(); config=Config()])

Plot summary panel of all fit-to-physio figures for MOC model
"""
function genfig_summary_fit_to_physio(
    params::MOCParams=MOCParams3(); 
    config=config, 
    color=colorschemes[:Dark2_8][1], 
    n_cf=9,
    level_step_probe=2.5,
    level_step_elicitor=2.5,
    kwargs...
)
    # Create figure and axes
    fig = Figure(; size=(1200, 400))
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    ax3 = Axis(fig[1, 3])
    ax4 = tradlogxax(fig[1, 4])

    # Run each figure and plot
    genfig_results_cas_threshold_vs_cf(params; fig=fig, ax=ax1, config=config, color=color, n_cf=n_cf, level_step_probe=level_step_probe, level_step_elicitor=level_step_elicitor, kwargs...)
    genfig_results_cas_ΔL_vs_cf(params; fig=fig, ax=ax2, config=config, color=color, n_cf=n_cf, level_step_probe=level_step_probe, level_step_elicitor=level_step_elicitor, kwargs...)
    genfig_results_cas_magnitude_vs_cf(params; fig=fig, ax=ax3, config=config, color=color, n_cf=n_cf, level_step_probe=level_step_probe, level_step_elicitor=level_step_elicitor, kwargs...)
    genfig_results_cas_sim_freq_sweeps(params; fig=fig, ax=ax4, config=config, kwargs...)

    # Add title
    toplab = "" 
    toplab *= "MOC efferent model simulation\n"
    toplab *= "width = $(params.moc_width)\n"
    toplab *= "offset = $(params.moc_offset)\n"
    toplab *= "beta = $(params.moc_beta)\n"
    Label(fig[0, 1]; tellwidth=false, text=toplab, halign=:left, fontsize=14.0)

    fig
end


"""
    genfig_results_cas_sim_freq_sweeps([params=MOCParams3(); config=Config()])

Plot iso-level tuning curves w/ and w/o a contralateral elicitor vs W&L 1989b Fig 6
"""
function genfig_results_cas_sim_freq_sweeps(
    params::MOCParams; 
    cf=2e3,
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogxax(fig[1, 1]),
    linewidth=2.0,
    markersize=5.0,
    level_elicitor=65.0,
    savefig=false,
)
    # Fetch model
    model = fetch_model(params, cf; stage="lsr")

    # Create simulations (iso-level TCs w/ probe at RLF midpoint)
    sims = map([-Inf, level_elicitor]) do lvl
        IsoLevelTC(
            model, 
            ProbePureToneElicitor3; 
            level_probe=lsr_midpoint(getcf(model); config=config), 
            level_elicitor=lvl
        )
    end

    # Add line at CF
    vlines!(ax, [getcf(model)/1e3]; color=:gray, linestyle=:dash, linewidth=linewidth/2)

    # Run and plot
    map(zip(sims,  [:black, :red])) do (sim, c)
        r = @memo config simulate(sim)
        scatter!(ax, Utilities.freq(sim) ./ 1e3, r; color=c, markersize=markersize)
        x̂, ŷ = smooth(log.(Utilities.freq(sim)), r, :loess; span=0.25)
        lines!(ax, exp.(x̂) ./ 1e3, ŷ; color=c, linewidth=linewidth)
    end

    # Add ticks and labels
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xlabel = "Elicitor frequency (kHz)"
    ax.ylabel = "Firing rate (sp/s)"

    # Adjust limits
    xlims!(ax, 0.1, 10.0)
    ylims!(ax, 0.0, 40.0)

    # Display and save
    display(fig)
    if savefig saveplot("fig_results_cas_sim_frequency_sweep_example.png", fig) end
    return fig, ax, sims
end

function genfig_results_cas_sim_freq_sweeps(
    params::Vector{<:MOCParams}; 
    cf=2e3,
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogxax(fig[1, 1]),
    linewidth=1.5,
    markersize=4.0,
    level_elicitor=65.0,
    colors=length(params) == 1 ? [:black] : vcat(:black, get.(Ref(colorschemes[:batlow]), LinRange(0.2, 1.0, length(params)))),
    savefig=false,
)
    # Loop over params
    map(zip(params, colors)) do (param, c)
        # Fetch model
        model = fetch_model(param, cf; stage="hsr")

        # Create simulations (iso-level TCs w/ probe at RLF midpoint)
        sim = IsoLevelTC(
            model, 
            ProbePureToneElicitor3; 
            level_probe=hsr_midpoint(getcf(model); config=config), 
            level_elicitor=level_elicitor
        )

        # Run and plot
        r = @memo config simulate(sim)
        scatter!(ax, Utilities.freq(sim) ./ 1e3, r; color=c, markersize=markersize)
        x̂, ŷ = smooth(log.(Utilities.freq(sim)), r, :loess; span=0.3)
        lines!(ax, exp.(x̂) ./ 1e3, ŷ; color=c, linewidth=linewidth)
    end

    # Add ticks and labels
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xlabel = "Elicitor frequency (kHz)"
    ax.ylabel = "Firing rate (sp/s)"

    # Adjust limits
    xlims!(ax, 0.1, 10.0)
    ylims!(ax, 75.0, 320.0)
    hlines!(ax, [85.0]; color=:gray, linestyle=:dash)
    text!(ax, [3.0], [87.0]; text="Spont", color=:gray)

    # Display and save
    display(fig)
    if savefig saveplot("fig_results_cas_sim_frequency_sweep_params.png", fig) end
    return fig, ax
end


"""
    genfig_results_cas_sim_cf_vs_bsf([params=[MOCParams3()]; config=Config()])
"""
function genfig_results_cas_sim_cf_vs_bsf(
    params::Vector{<:MOCParams}=[MOCParams3()]; 
    cf=LogRange(0.5e3, 16e3, 15),
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogdoubleax(fig[1, 1]),
    markersize=5.0,
    level_elicitor=70.0,
    colors=length(params) == 1 ? [:black] : vcat(:black, get.(Ref(colorschemes[:batlow]), LinRange(0.0, 1.0, length(params)))),
    savefig=false,
)
    # Draw abline
    x = LogRange(0.1, 30.0, 1000)
    lines!(ax, x, x; color=:gray, linewidth=2.0)

    # Select elicitor freqs and annotate (very faintly, just for debug) what's going on
    freqs_elicitor = LogRange(0.18e3, 8e3, 41)
    #hlines!(ax, freqs_elicitor ./ 1e3; color=(:gray, 0.05), linewidth=0.5)

    # Loop over params
    map(zip(params, colors)) do (param, c)
        bsfs = map(cf) do _cf
            # Fetch model
            model = fetch_model(param, _cf)

            # Create simulations (iso-level TCs w/ probe at RLF midpoint)
            sim = IsoLevelTC(
                model, 
                ProbePureToneElicitor3; 
                level_probe=hsr_midpoint(getcf(model); config=config), 
                level_elicitor=level_elicitor, 
                freqs_elicitor=freqs_elicitor,
            )

            # Run and plot
            r = @memo config simulate(sim)
            idxmin = argmin(r)
            return Utilities.freq(sim)[idxmin]
        end
        scatter!(ax, cf ./ 1e3, bsfs ./ 1e3; color=c, markersize=markersize)
    end

    # Add ticks and labels
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "BSF (kHz)"

    # Adjust limits
    xlims!(ax, 0.1, 10.0)
    ylims!(ax, 0.1, 10.0)

    # Display and save
    display(fig)
    if savefig saveplot("fig_results_cas_sim_cf_vs_bsf.png", fig) end
    return fig, ax
end

"""
    genfig_results_cas_sim_cf_vs_bsf([params=[MOCParams3()]; config=Config()])
"""
function genfig_results_cas_sim_cf_vs_bsf_analyze(
    params::Vector{<:MOCParams}=[MOCParams3()]; 
    cf=[1e3, 2e3, 3e3, 4e3],
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=tradlogdoubleax(fig[1, 1]),
    markersize=5.0,
    level_elicitor=70.0,
    colors=length(params) == 1 ? [:black] : vcat(:black, get.(Ref(colorschemes[:batlow]), LinRange(0.0, 1.0, length(params)))),
)
    # Draw abline
    x = LogRange(0.1, 30.0, 1000)
    lines!(ax, x, x; color=:gray, linewidth=2.0)

    # Select elicitor freqs and annotate (very faintly, just for debug) what's going on
    freqs_elicitor = LogRange(0.5e3, 8e3, 21)
    hlines!(ax, freqs_elicitor ./ 1e3; color=(:gray, 0.05), linewidth=0.5)

    # Loop over params
    map(zip(params, colors)) do (param, c)
        bsfs = map(cf) do _cf
            # Fetch model
            model = fetch_model(param, _cf; species="human")

            # Create simulations (iso-level TCs w/ probe at RLF midpoint)
            sim = IsoLevelTC(
                model, 
                ProbePureToneElicitor3; 
                level_probe=hsr_midpoint(getcf(model); config=config), 
                level_elicitor=level_elicitor, 
                freqs_elicitor=freqs_elicitor,
            )

            # Run and plot
            r = @memo config simulate(sim)
            idxmin = argmin(r)
            return Utilities.freq(sim)[idxmin]
        end
        scatter!(ax, cf ./ 1e3, bsfs ./ 1e3; color=c, markersize=markersize)
    end

    # Add ticks and labels
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "BSF (kHz)"

    # Adjust limits
    xlims!(ax, 0.1, 10.0)
    ylims!(ax, 0.1, 10.0)

    # Display and save
    display(fig)
    return fig, ax
end


"""
    genfig_results_cas_band_widening([params=[MOCParams3()]; config=Config()])
"""
function genfig_results_cas_band_widening(
    params::Vector{<:MOCParams}=[MOCParams3()]; 
    cf=2e3,
    config=Config(),
    size=(2.0*100, 1.6*100),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    markersize=5.0,
    level_elicitor=20.0,
    colors=length(params) == 1 ? [:black] : vcat(:black, get.(Ref(colorschemes[:batlow]), LinRange(0.2, 1.0, length(params)))),
    savefig=false,
)
    # Create stimuli (probe tone with noise elicitor)
    probe = PaddedStimulus(
        PureTone(; freq=cf, dur=0.1, level=hsr_midpoint(cf; config=config), fs=100e3),
        0.9,
        0.0,
        100e3
    )
    bandwidths = 1/3 .* (1:1:8)
    stimuli = map(bandwidths) do bandwidth
        elicitor = GaussianNoiseSL(; 
            freq_low=cf * 2^(-bandwidth/2), 
            freq_high=cf * 2^(bandwidth/2), 
            level=level_elicitor, 
            dur=1.0,
            fs=100e3,
        )
        BinauralStimulus(probe, elicitor)
    end

    # Create control stimulus (probe tone with tone elicitor, level matching 1/3 octave noise overall level)
    probe = PaddedStimulus(
        PureTone(; freq=cf, dur=0.1, level=hsr_midpoint(cf; config=config), fs=100e3),
        0.9,
        0.0,
        100e3
    )
    elicitor = PureTone(; 
        freq=cf,
        level=dbspl(synthesize(contra(stimuli[1]))), 
        dur=1.0,
        fs=100e3,
    )
    stimulus_control = BinauralStimulus(probe, elicitor)

    # Fetch model and plot how rate changes vs elicitor bandwidth
    map(zip(params, colors)) do (param, c)
        # Fetch model
        model = fetch_model(param, cf)

        # Run all stims, averaging rates from 0.9-1.0 s matching other sims
        r = map(stimuli) do stim
            resp = @memo config compute(model, stim)
            mean(resp[samples(0.9, 1.0, samprate(model))])
        end
        scatter!(ax, bandwidths, r; color=c, markersize=markersize)
        x̂, ŷ = smooth(bandwidths, r, quadratic)
        lines!(ax, x̂, ŷ; color=c, linewidth=2.0)
    end

    # Plot control result
    map(zip(params, colors)) do (param, c)
        # Fetch model
        model = fetch_model(param, cf)

        # Run control stimulus, averaging rates from 0.9-1.0 s matching other sims
        resp = @memo config compute(model, stimulus_control)
        r = mean(resp[samples(0.9, 1.0, samprate(model))])
        scatter!(ax, [-0.5], [r]; color=c, markersize=markersize)
    end

    # Add line for spont
    hlines!(ax, [85.0]; color=:gray, linestyle=:dash, linewidth=1.0)
    text!(ax, [2.0], [87.0]; text="Spont", color=:gray)

    # Limits
    ylims!(ax, 80.0, 220.0)
    xlims!(ax, -1.0, 3.0)

    # Labels
    ax.xlabel = "Noise bandwidth (oct)"
    ax.ylabel = "Firing rate (sp/s)"
    ax.xminorticksvisible = true
    ax.xticks = ([0.0, 1.0, 2.0, 3.0], ["", "1", "2", "3"])
    ax.xminorticks = IntervalsBetween(3)

    # Display and save
    display(fig)
    if savefig saveplot("fig_results_cas_band_widening.png", fig) end
    return fig, ax
end

# Function to plot RLF shift results in elegant format
function genfig_results_summary_rlf(
    params::MOCParams=MOCParams3(); 
    config=Config(),
    cf=2e3,
    levels=-10.0:2.5:80.0,
    levels_elicitor=vcat(-Inf, 0.0:20.0:80.0),
    durs_elicitor=0.02 .* 2 .^ (0:2:6),
    gaps_elicitor=[0.0, 0.1, 0.2, 0.4],
    fig=Figure(; size=(8.0*100, 2.5*100)),
    colors=cs_level.(levels_elicitor),
    inset_xpos=1.0,
    inset_ypos=0.45,
)
    # Select model for simulations
    model = fetch_model(params, cf)

    # Create axes and inset axes for level plot
    ax_level = Axis(
        fig[1, 1];
        xlabel="Probe level (dB SPL)",
        ylabel="Firing rate (sp/s)",
        title="Effect of elicitor level"
    )
    xlims!(ax_level, -5.0, 55.0)
    ylims!(ax_level, 0.0, 360.0)

    ax_level_inset = Axis(
        fig[1, 1],
        width=Relative(0.28),
        height=Relative(0.28),
        halign=inset_xpos,
        valign=inset_ypos,
        xlabel="E. level (dB SPL)",
        ylabel="ΔL (dB SPL)",
    )
    xlims!(ax_level_inset, -10.0, 90.0)
    ylims!(ax_level_inset, -2.0, 15.0)

    # Create axes and inset axes for elicitor duration plot
    ax_dur = Axis(
        fig[1, 2];
        xlabel="Probe level (dB SPL)",
        ylabel="Firing rate (sp/s)",
        title="Effect of elicitor duration"
    )
    xlims!(ax_dur, -5.0, 55.0)
    ylims!(ax_dur, 0.0, 360.0)

    ax_dur_inset = Axis(
        fig[1, 2],
        width=Relative(0.28),
        height=Relative(0.28),
        halign=inset_xpos,
        valign=inset_ypos,
        xlabel="E. duration (ms)",
        ylabel="ΔL (dB SPL)",
        xscale=log10,
        xticks=durs_elicitor .* 1e3,
        xticklabelrotation=π/4,
    )
    xlims!(ax_dur_inset, minimum(durs_elicitor)/2*1e3, maximum(durs_elicitor)*2*1e3)
    ylims!(ax_dur_inset, -2.0, 15.0)

    # Create axes and inset axes for elicitor-probe gap
    ax_gap = Axis(
        fig[1, 3];
        xlabel="Probe level (dB SPL)",
        ylabel="Firing rate (sp/s)",
        title="Effect of elicitor-probe gap"
    )
    xlims!(ax_gap, -5.0, 55.0)
    ylims!(ax_gap, 0.0, 360.0)

    ax_gap_inset = Axis(
        fig[1, 3],
        width=Relative(0.28),
        height=Relative(0.28),
        halign=inset_xpos,
        valign=inset_ypos,
        xlabel="E.-P. gap (ms)",
        ylabel="ΔL (dB SPL)",
        xscale=log10,
        xticks=durs_elicitor .* 1e3,
    )
    xlims!(ax_dur_inset, minimum(durs_elicitor)/2*1e3, maximum(durs_elicitor)*2*1e3)
    ylims!(ax_gap_inset, -2.0, 15.0)

    # Map over elicitor levels and fill in elicitor-level plot
    map(zip(levels_elicitor, colors)) do (elvl, c)
        sim = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=elvl)
        r = @memo config simulate(sim)
        lines!(ax_level, level(sim), r; color=c, linewidth=2.5)
    end
    ΔL = map(levels_elicitor[2:end]) do elvl
         sim_np = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=-Inf)
         sim_wp = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=elvl)
         return threshold_shift(sim_np, sim_wp; config=config)
    end
    lines!(ax_level_inset, levels_elicitor[2:end], ΔL; color=:gray)
    scatter!(ax_level_inset, levels_elicitor[2:end], ΔL; color=colors[2:end])

    # Map over elicitor durations and fill in elicitor-dur plot
    colors = colorschemes[:Dark2_8]
    map(zip(durs_elicitor, colors)) do (dur, c)
        sim = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=60.0, dur_delay=-dur, dur_suppressor=0.0, dur_wait=maximum(durs_elicitor)-dur)
        r = @memo config simulate(sim)
        lines!(ax_dur, level(sim), r; color=c, linewidth=2.5)
    end
    ΔL = map(durs_elicitor) do dur
         sim_np = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=-Inf, dur_delay=-dur, dur_suppressor=0.0, dur_wait=maximum(durs_elicitor)-dur)
         sim_wp = RLF(model, ProbePureToneElicitor3; levels=levels, level_elicitor=60.0, dur_delay=-dur, dur_suppressor=0.0, dur_wait=maximum(durs_elicitor)-dur)
         return threshold_shift(sim_np, sim_wp; config=config)
    end
    lines!(ax_dur_inset, durs_elicitor .* 1e3, ΔL; color=:gray)
    scatter!(ax_dur_inset, durs_elicitor .* 1e3, ΔL; color=colors[1:length(durs_elicitor)])

    # # Map over elicitor durations and fill in elicitor-dur plot
    # colors = colorschemes[:Dark2_8]
    # map(zip(gaps_elicitor, colors)) do (gap, c)
    #     sim = RLF(model, ProbePureToneElicitor; levels=levels, level_elicitor=60.0, dur??=gap)
    #     r = @memo config simulate(sim)
    #     lines!(ax_dur, level(sim), r; color=c, linewidth=2.5)
    # end
    # ΔL = map(durs_elicitor) do dur
    #      sim_np = RLF(model, ProbePureToneElicitor; levels=levels, level_elicitor=-Inf, dur??=gap)
    #      sim_wp = RLF(model, ProbePureToneElicitor; levels=levels, level_elicitor=60.0, dur??=gap)
    #      return threshold_shift(sim_np, sim_wp; config=config)
    # end
    # lines!(ax_dur_inset, gaps_elicitor .* 1e3, ΔL; color=:gray)
    # scatter!(ax_dur_inset, gaps_elicitor .* 1e3, ΔL; color=colors[1:length(gaps_elicitor)])

    # Plot and save
    fig
end

# ##########################################################################################
# MORE SUPPLEMENTARY FIGURES
# - genfig_maximum_attenuation()
function genfig_maximum_attenuation(params; cfs=LogRange(1e3, 16e3, 15), config=Config(), savefig=true)
    # First, map over all CFs and simulate gain reduction elicited by a 100 dB SPL gaussian noise
    ΔLs_100dB = pmap(cfs) do cf
        model = fetch_model(params, cf; stage="gainpostmix")
        stim = GaussianNoise(; dur=1.0, level=100.0)
        cohc_to_ΔL(minimum(compute(model.model_contra, stim)); cf=cf, config=config)
    end
    println("Max attn at 100 dB SPL elicitor: $(maximum(ΔLs_100dB)) dB")

    # Second, compute what the ΔL should be if a spike rate of 100 sp/s was really delivered by
    # the MOC system, before gain-mixing
    ŷ = moc_nonlinearity.(100.0, Ref(params), cfs)
    ΔL_100sp = pmap(zip(cfs, ŷ)) do (x, y)
        cohc_to_ΔL(y; cf=x, config=config)
    end
    println("Max attn at 100 sp/s: $(maximum(ΔL_100sp)) dB")

    # Create two-panel figure to show everything
    fig = Figure(; size=(5*100, 2.5*100))
    ax1 = tradlogxax(fig[1, 1], xlabel="CF (kHz)", ylabel="ΔL (dB)", xscale=log10)
    ax2 = tradlogxax(fig[1, 2], xlabel="CF (kHz)", ylabel="ΔL (dB)", xscale=log10)
    ax1.xticks = [0.1, 1, 10]
    ax2.xticks = [0.1, 1, 10]
    scatter!(ax1, cfs ./ 1e3, ΔLs_100dB; color=:black, markersize=10)
    scatter!(ax2, cfs ./ 1e3, ΔL_100sp; color=:black, markersize=10)
    xlims!(ax1, 0.1, 30.0)
    ylims!(ax1, 0.0, 55.0)
    xlims!(ax2, 0.1, 30.0)
    ylims!(ax2, 0.0, 55.0)

    display(fig)
    if savefig
        saveplot("fig_supp_maximum_attenuation.png", fig) 
    end
 
    fig
end

function genfig_followup_elicitor_eps(params; config=Config(), level_elicitor=65.0)
    # Fetch model 
    model = fetch_model(params, 2e3)

    # Create simulations (iso-level TCs w/ probe at RLF midpoint)
    sim = IsoLevelTC(
        model, 
        ProbePureToneElicitor3; 
        level_probe=hsr_midpoint(getcf(model); config=config), 
        level_elicitor=level_elicitor,
    )

    # Extract all the contralateral elicitors from the simulation
    elicitors = map(sim.probes) do stim
        contra(stim)
    end

    # Plot elicitor responses near versus off-CF
    f_elicitor = Utilities.freq.(elicitors)
    idx_min = argmin(abs.(f_elicitor .- getcf(model)))
    
    # Create figure window
    fig = Figure(; size=(1000, 300))
    ax = Axis(fig[1, 1]; xscale=log10)
    ax_premix = Axis(fig[1, 2]; xscale=log10)
    ax_postmix = Axis(fig[1, 3]; xscale=log10)
    axs = [ax, ax_premix, ax_postmix]
    map(axs) do ax
        xlims!(ax, (extrema(getcf(contra(model))) ./ 1e3)...)
        ax.xticks = [0.5e3, 1e3, 2e3, 4e3, 8e3] ./ 1e3
        ax.xlabel = "Contralateral CF (kHz)"
    end
    ylims!(ax, 0.0, 80.0)
    ylims!(ax_premix, 0.0, 1.0)
    ylims!(ax_postmix, 0.0, 1.0)
    ax.ylabel = "LSR average firing rate (sp/s)"
    ax_premix.ylabel = "Premix γ [0, 1]"
    ax_postmix.ylabel = "Postmix γ [0, 1]"
    fig

    function normpdf(x, mu, sigma)
        var = sigma ^ 2.0
        normfac = 1/sqrt(2π * var);
        return normfac * exp(-1 * (x-mu)^2.0 / (2*var));
    end

    # Note about what is happening to compute postmix gain:
    # for (int c=0; c < n_chan; c++) {
    #     double weight_total = 0.0;
    #     for (int subc = 0; subc < n_chan; subc++) {
    #         gainpostmix[c][n] = gainpostmix[c][n] * pow(gain[subc][n], moc_wb_weights[c][subc]);
    #         weight_total = weight_total + moc_wb_weights[c][subc];
    #     }
    #     gainpostmix[c][n] = pow(gainpostmix[c][n], 1 / weight_total);
    # }
    weight_curve = normpdf.(log2(getcf(model)), log2.(getcf(contra(model))), 0.9)

    # Simulate near-CF response
    colors = [:black, :red, :pink]
    idxs = [idx_min, idx_min-1, idx_min-2]
    map(zip(idxs, colors)) do (idx, c)
        # Plot neural responses to elicitor in LSR contra pathway
        r_full = @memo config Utilities._compute(model.model_contra, elicitors[idx])
        r = map(mean, r_full["lsr"])
        lines!(ax, getcf(contra(model)) ./ 1e3, r; color=c, linewidth=2.0)

        # Add text label at top of each curve indicating sum of rates
        total_rate = sum(r)
        text!(ax, [f_elicitor[idx] / 1e3], [50.0]; text="$(round(total_rate))", color=c)

        # Plot premix gain reduction arising from each channel in response to elicitor
        r = map(mean, r_full["gain"])
        lines!(ax_premix, getcf(contra(model)) ./ 1e3, r; color=c, linewidth=2.0)
        vlines!(ax_premix, [2.0]; color=:gray, linestyle=:dash)
        lines!(ax_premix, getcf(contra(model)) ./ 1e3, weight_curve; color=:gray, linestyle=:dash)
        
        # Add text label at bottom of each curve indicating sum of gain reduction
        total_gain_reduction = sum(-1.0 .* (r .- 1.0))
        text!(ax_premix, [f_elicitor[idx] / 1e3], [0.3]; text="$(round(total_gain_reduction, digits=2))", color=c)

        # Plot postmix gain reduction in response to elicitor
        r = map(mean, r_full["gainpostmix"])
        lines!(ax_postmix, getcf(contra(model)) ./ 1e3, r; color=c, linewidth=2.0)
        vlines!(ax_postmix, [2.0]; color=:gray, linestyle=:dash)
        lines!(ax_premix, getcf(contra(model)) ./ 1e3, weight_curve; color=:gray, linestyle=:dash)
    end

    fig
end