# WarrenLiberman1989.jl
#
# Source code for simulating experiments from Warren and Liberman (1989):
# Warren, E. H., & Liberman, M. C. (1989a). Effects of contralateral sound on auditory-nerve 
# responses. I. Contributions of cochlear efferents. Hearing Research, 37(2), 89–104. 
# https://doi.org/10.1016/0378-5955(89)90032-4
# Warren, E. H., & Liberman, M. C. (1989b). Effects of contralateral sound on auditory-nerve 
# responses. II. Dependence on stimulus variables. Hearing Research, 37(2), 105–121. 
# https://doi.org/10.1016/0378-5955(89)90033-6

export WarrenTimeCourse, WarrenTuningCurve, analyze_full, annotate, freq_probe, level_probe

include(joinpath("data.jl"))

"""
    WarrenTimeCourse <: Simulation

Simulation of time course of ipsilateral response during CAS.

Provides the following methods: 
- `WarrenTimeCourse(cf)`: Create time-course simulation at specified CF
- `viz(sim)`: visualize main and supplemental outputs with no caching
- `neurogram(sim)`: visualize contralateral population response with no caching
- `viz(sims::Vector[, ::Config])`: visualize multiple time courses together with caching
- `idxswin(sim)`: return indices for two analysis windows
- `analyze(sim)`: return rate change due to CAS in percent, without caching
- `analyze(sim, ::Config)`: return rate change due to CAS in percent, with caching
- `threshold(sims::Vector[, ::Config])`: return level at which rate change surpasses criterion
"""
@with_kw struct WarrenTimeCourse{A, B} <: Simulation where {A <: AbstractStimulus, B <: AbstractStimulus}
    stim_ipsi::A    
    stim_contra::B
    stim::BinauralStimulus=BinauralStimulus(stim_ipsi, stim_contra)
    model::GFC2024_FakeBin
end

# Default external constructor
function WarrenTimeCourse(
    cf; 
    dur=2.0, 
    dur_sup=1.0,
    delay_sup=1.0,
    level=NaN, 
    level_sup_re_probe=25.0,
    sup=GaussianNoise,
    freq_sup=cf,
    kwargs...
)
    # Create model 
    model = GFC2024_FakeBin(cf; cf_ipsi=cf, kwargs...)

    # If level is NaN, we automatically determine the level for probe based on identifying
    # the midpoint of the ipsilateral rate-level function
    if isnan(level)
        level = midpoint(RLF(ipsi(model)))
    end
    stim_ipsi = PureTone(; freq=cf, dur=dur, level=level)

    # Create contralateral stimulus
    if sup == GaussianNoise
        temp = GaussianNoise(; dur=dur_sup, level=level+level_sup_re_probe)
    elseif sup == PureTone
        temp = PureTone(; freq=freq_sup, dur=dur_sup, level=level+level_sup_re_probe)
    else
        error("Suppresor stimulus type not recognized!")
    end
    stim_contra = MixedStimulus(
        Silence(; dur=dur),
        temp;
        d=delay_sup,
    )

    # Combine and bundle into WarrenTimeCourse
    stim = BinauralStimulus(stim_ipsi, stim_contra)
    WarrenTimeCourse(; stim_ipsi=stim_ipsi, stim_contra=stim_contra, stim=stim, model=model)
end

# `freq` method (points to contralateral frequency)
Utilities.freq(sim::WarrenTimeCourse) = Utilities.freq(contra(stimobj(sim)).b)
freq_probe(sim::WarrenTimeCourse) = Utilities.freq(ipsi(stimobj(sim)))

# `level` method (points to contralateral level)
Utilities.level(sim::WarrenTimeCourse) = Utilities.level(contra(stimobj(sim)).b)
level_probe(sim::WarrenTimeCourse) = Utilities.level(ipsi(stimobj(sim)))

# `simulate` method
function Utilities.simulate(sim::WarrenTimeCourse)
    compute(sim.model, sim.stim)
end

# `idxswin` method for the simulation itself, returning idxs for two windows:
#   1) Post probe onset, after rapid onset adaptation but before suppressor begins
#   2) Post suppressor onset, after ~100 ms of buildup is allowed to take place
function Utilities.idxswin(
    sim::WarrenTimeCourse; 
    dur_skip_onset=min(dur(ipsi(stimobj(sim)))/4, 0.4),
    dur_wait_moc=min(dur(contra(stimobj(sim)).b)/4, 0.4),
)
    stim = sim.stim
    idx_onset_probe = sampleat(onset(ipsi(stim)), samprate(sim))
    idx_start = idx_onset_probe + samples(dur_skip_onset, samprate(sim))
    idxs_probe = idx_start:samples(contra(stim).d, samprate(sim))
    idx_onset_sup = sampleat(contra(stim).d, samprate(sim))
    idx_start = idx_onset_sup + samples(dur_wait_moc, samprate(sim))
    idxs_sup = idx_start:(samples(contra(stim).d, samprate(sim)) + length(contra(stim).b))
    return idxs_probe, idxs_sup
end

# `analyze` method for single time course, returning change in rate in percent
Utilities.analyze(sim::WarrenTimeCourse; kwargs...) = analyze(sim, simulate(sim); kwargs...)
Utilities.analyze(sim::WarrenTimeCourse, config::Config; kwargs...) = analyze(sim, @memo config simulate(sim); kwargs...)
function Utilities.analyze(sim::WarrenTimeCourse, r::Vector{Float64}; kwargs...)
    # Determine rates in two windwos of interest and compute difference (in %)
    idxs_probe, idxs_sup = idxswin(sim; kwargs...)
    μ_probe = mean(r[idxs_probe])
    μ_sup = mean(r[idxs_sup])
    return 100.0 * (μ_sup - μ_probe)/μ_probe
end

# `analyze` method for single time course, returning change in rate and time constants
analyze_full(sim::WarrenTimeCourse; kwargs...) = analyze_full(sim, simulate(sim); kwargs...)
analyze_full(sim::WarrenTimeCourse, config::Config; kwargs...) = analyze_full(sim, @memo config simulate(sim); kwargs...)
function analyze_full(
    sim::WarrenTimeCourse, 
    r::Vector{Float64}; 
    dur_window=4e-3, 
    dur_skip_onset=min(dur(ipsi(stimobj(sim)))/4, 0.4),
    dur_wait_moc=min(dur(contra(stimobj(sim)).b)/4, 0.4),
    kwargs...
)
    # First, calculate Δ
    Δ = analyze(sim, r; dur_skip_onset=dur_skip_onset, dur_wait_moc=dur_wait_moc, kwargs...)

    # Next, also characterize time courses
    # In Warren and Liberman (1989), they specify the time courses in terms of time needed 
    # to achieve 90% of maximal suppression (for onset) and time to achieve 90% of recovery
    # First, calculate smoothed time points and rates (eliminate phase locking)
    t_smooth, r_smooth = slidwin(r; dur_window=dur_window, kwargs...)  # smooth rate response

    # Next, calculate rates within relevant time windows (from smoothed rates) 
    idxs = idxswin(sim; dur_skip_onset=dur_skip_onset, dur_wait_moc=dur_wait_moc)               
    times_probe = [timeat(idxs[1][1], samprate(sim)), timeat(idxs[1][end], samprate(sim))]
    times_sup = [timeat(idxs[2][1], samprate(sim)), timeat(idxs[2][end], samprate(sim))]
    μ_probe = mean(r_smooth[(t_smooth .>= times_probe[1]) .& (t_smooth .<= times_probe[2])])  # rate during probe-only window
    μ_sup = mean(r_smooth[(t_smooth .>= times_sup[1]) .& (t_smooth .<= times_sup[2])])

    # Finally, identify 1) time needed to achieve 90% of maximal suppression...
    t_onset_sup = onset(contra(stimobj(sim)))
    μ_90pct_sup = μ_probe + 0.9*Δ/100 * μ_probe
    rate_achieved = r_smooth .< μ_90pct_sup
    rate_achieved[t_smooth .< t_onset_sup] .= 0
    idxs = findall(rate_achieved)
    τ_onset = isempty(idxs) ? NaN : round(t_smooth[first(idxs)] - t_onset_sup; digits=3)  # round to ms precision

    # ... and 2) time needed to recover 90%
    t_offset_sup = offset(contra(stimobj(sim)))
    μ_90pct_rec = μ_probe + 0.1*Δ/100 * μ_probe
    rate_achieved = r_smooth .> μ_90pct_rec
    rate_achieved[t_smooth .< t_offset_sup] .= 0
    idxs = findall(rate_achieved)
    τ_offset = isempty(idxs) ? NaN : round(t_smooth[first(idxs)] - t_offset_sup; digits=3)  # round to ms precision

    # Return all
    return Δ, τ_onset, τ_offset, μ_probe, μ_sup
end

# `threshold` for series of time courses, determines a threshold for suppressor level vs
# suppression percentage curve
function Utilities.threshold(sims::Vector{<:WarrenTimeCourse}, config=Config(); criterion=-5.0)
    # Get levels from all simulations (level of contra stimulus, dB SPL)
    levels = level.(sims)

    # Assert that orders are in level
    @assert issorted(levels)

    # Get suppression amounts from each curve
    Δs = analyze.(sims, Ref(config))

    # Subtract out base amount drift (we need a level at least below 20 dB SPL to do this)
    Δ_baseline = mean(Δs[levels .<= 20.0])
    Δs .-= Δ_baseline

    # Determine threshold as first point to cross criterion
    idxs_thr = findall(x -> x < criterion, Δs)
    return isempty(idxs_thr) ? NaN : levels[first(idxs_thr)]
end

# `viz` methods
function Utilities.viz(sim::WarrenTimeCourse)
    viz(sim.model, sim.stim)
end

Utilities.viz(sim::WarrenTimeCourse, config::Config; kwargs...) = viz([sim], config; include_legend=false, kwargs...)

function Utilities.viz(
    sims::Vector{<:WarrenTimeCourse}, 
    config=Config();
    label="MOC weight",
    include_legend=true,
    include_inline_legend=false,
    labeler=x -> string(contra(x.model).moc_weight_wdr),
    colors=length(sims) == 1 ? [:black] : get(colorschemes[:Dark2_8], LinRange(0.0, 1.0, length(sims))),
    plot_windows=true,
    plot_time_constants=false,
    size=(800, 450),
    ylims=(0, 350.0),
    linewidth=2.0,
    kwargs...
)
    # Get curves from simulations
    curves = pmap(enumerate(sims)) do (idx, sim)
        @memo config simulate(sim)
    end

    # Create figure and axes
    fig = Figure(; size=size)
    ax = Axis(fig[1, 1])

    # Add markers to indicate time windows
    if plot_windows
        idxs_probe, idxs_sup = idxswin(sims[1]; kwargs...)
        t_probe = collect(timeat.(extrema(idxs_probe), samprate(sims[1])))
        t_sup = collect(timeat.(extrema(idxs_sup), samprate(sims[1])))
        vspan!(ax, collect.(zip(t_probe, t_sup))...; color=(:gray, 0.25))
        text!(ax, [t_probe[1]], [ylims[2] - 5.0]; text="Probe\nwindow", rotation=π/2, align=(:right, :top), color=:white, font=:bold)
        text!(ax, [t_sup[1]], [ylims[2] - 5.0]; text="Elicitor\nwindow", rotation=π/2, align=(:right, :top), color=:white, font=:bold)
    end

    # Plot each curve as line
    lns = map(zip(sims, curves, colors)) do (sim, curve, c)
        Δ, τ_onset, τ_offset, μ_onset, μ_offset = analyze_full(sim, config)
        lines!(ax, slidwin(curve; dur_window=4e-3)...; label="$(labeler(sim))", color=c, linewidth=linewidth)
    end

    # Create labels for legend, and optionally annotate plot with markers for time course recovery
    labels = map(zip(sims, colors)) do (sim, c)
        # Analyze results 
        Δ, τ_onset, τ_offset, μ_onset, μ_offset = analyze_full(sim, config)

        # Annotate
        if plot_time_constants
            annotate(ax, sim, config; color=color)
        end

        # Create and return label
        "$(labeler(sim))" # → $(round(Δ; digits=2))%"
    end

    # Add legend 
    if include_legend
        Legend(fig[1, 2], lns, labels, label; fontsize=10.0)
    end
    if include_inline_legend
        axislegend(; position=:rt, orientation=:horizontal, nbanks=2)
    end

    # Adjust axes and limits
    ax.xlabel = "Time (s)"
    ax.ylabel = "Firing rate (sp/s)"
    ylims!(ax, ylims...)
    xticks = 0.0:0.5:dur(sims[1].stim)
    xticklabels = string.(xticks)
    xticklabels[2:2:end] .= ""
    ax.xticks = (xticks, xticklabels)

    # Return figure
    fig
end

# `annotate` method defined on axis
function annotate(
    ax::Axis, 
    sim::WarrenTimeCourse, 
    config::Config;
    linestyle=:dash,
    color=:gray,
    )
    # Perform analysis
    Δ, τ_onset, τ_offset, μ_onset, μ_offset = analyze_full(sim, config)

    # Plot lines to indicate onset time constant
    lines!(ax, [onset(contra(stimobj(sim))), onset(contra(stimobj(sim))) + τ_onset], [μ_onset, μ_onset]; color=color, linestyle=linestyle)
    lines!(ax, [onset(contra(stimobj(sim))) + τ_onset, onset(contra(stimobj(sim))) + τ_onset], [μ_onset, 0.0]; color=color, linestyle=linestyle)

    # Plot lines to indicate offset time constant
    lines!(ax, [offset(contra(stimobj(sim))), offset(contra(stimobj(sim))) + τ_offset], [μ_offset, μ_offset]; color=color, linestyle=linestyle)
    lines!(ax, [offset(contra(stimobj(sim))) + τ_offset, offset(contra(stimobj(sim))) + τ_offset], [μ_offset, 0.0]; color=color, linestyle=linestyle)

    # Add text
    if !isnan(τ_onset)
        text!(ax, [onset(contra(stimobj(sim))) + τ_onset], [μ_onset * 1.1]; text="τ = $(Int(round(τ_onset*1e3))) ms", color=color, fontsize=10.0)
    end
    if !isnan(τ_offset)
        text!(ax, [offset(contra(stimobj(sim))) + τ_offset + 0.05], [μ_offset * 0.4]; text="τ = $(Int(round(τ_offset*1e3))) ms", color=color, fontsize=10.0)
    end
end

# `neurogram` method
function Utilities.neurogram(sim::WarrenTimeCourse)
    neurogram(sim.model, sim.stim)
end

"""
    WarrenTuningCurve <: Simulation

Simulation of frequency profile of CAS with varying suppressor frequency

Provides the following methods: 
- `WarrenTuningCurve(cf)`: Create frequency-tuning simulation at specified CF
- `simulate(sim)`: simulate without caching by evaluating `simulate` on each constituent sim
- `simulate(sim, ::Config)`: simulate with caching
- `freq(sim)`: return frequency of contralateral suppressor
- `viz(sim[, ::Config])`: visualize profile with caching
"""
@with_kw struct WarrenTuningCurve <: Simulation
    model::GFC2024_FakeBin
    sims::Vector{WarrenTimeCourse{PureTone, MixedStimulus{Silence, PureTone}}}
end

# Default external constructor
function WarrenTuningCurve(
    cf; 
    dur=3.0, 
    dur_sup=1.0,
    delay_sup=1.0,
    level=NaN, 
    level_sup_re_probe=20.0,
    freqs_sup=LogRange(0.1e3, 20e3, 15),
    cf_low=0.2e3,
    cf_high=10e3,
    n_cf=41,
    kwargs...
)
    # Create model 
    model = GFC2024_FakeBin(LogRange(cf_low, cf_high, n_cf); cf_ipsi=cf, kwargs...)

    # Create ipsilateral stimulus
    # If level is NaN, we automatically determine the level for probe based on identifying
    # the midpoint of the ipsilateral rate-level function
    if isnan(level)
        level = midpoint(RLF(ipsi(model)))
    end
    stim_ipsi = PureTone(; freq=cf, dur=dur, level=level)

    # Create contralateral stimuli and bundle everything into simulation objects
    sims = map(freqs_sup) do freq_sup
        temp = PureTone(; freq=freq_sup, dur=dur_sup, level=level+level_sup_re_probe)
        stim_contra = MixedStimulus(
            Silence(; dur=dur),
            temp;
            d=delay_sup,
        )
        stim = BinauralStimulus(stim_ipsi, stim_contra)
        WarrenTimeCourse(; stim_ipsi=stim_ipsi, stim_contra=stim_contra, stim=stim, model=model)
    end

    # Create
    WarrenTuningCurve(model, sims)
end

# Small generic methods
Utilities.freq(sim::WarrenTuningCurve) = getfield.(getfield.(contra.(getfield.(sim.sims, :stim)), :b), :freq)

# `simulate` methods: one for simulating each curve fresh, another for using cached copies internally
function Utilities.simulate(sim::WarrenTuningCurve; kwargs...)
    # Simulate each time course and get Δ
    @showprogress pmap(sim.sims) do _sim
        analyze(_sim; kwargs...)
    end
end

function Utilities.simulate(sim::WarrenTuningCurve, config::Config; kwargs...)
    # Simulate each time course and get Δ
    @showprogress pmap(enumerate(sim.sims)) do (idx, _sim)
        analyze(_sim, config; kwargs...)
    end
end

# `analyze` methods
function Utilities.analyze(sim::WarrenTuningCurve, config::Config; kwargs...)
    r = simulate(sim, config; kwargs...)
    Utilities.analyze(sim, r)
end

function Utilities.analyze(sim::WarrenTuningCurve, r::Vector{Float64})
    # Identify bsf (best suppression frequency) and maximum suppression
    bsf(sim, r), minimum(r)
end

bsf(sim::WarrenTuningCurve, r::Vector{Float64}) = Utilities.freq(sim)[argmin(r)]

function Utilities.analyze(sims::Vector{WarrenTuningCurve}, config::Config; kwargs...)
    map(sims) do sim
        Utilities.analyze(sim, config; kwargs...)
    end
end

# `viz` method
function Utilities.viz(sim::WarrenTuningCurve, config::Config, fig::Figure=Figure(); kwargs...)
    ax = Axis(
        fig[1, 1]; 
        xscale=log10, 
        xminorticksvisible=true, 
        xminorticks=IntervalsBetween(9),
        xticks=[0.1, 1.0, 10.0],
    )
    viz(sim, config, fig, ax; kwargs...)
end

function Utilities.viz(
    sims::Vector{<:WarrenTuningCurve}, 
    config::Config, 
    fig::Figure=Figure(); 
    colors=colorschemes[:Dark2_8],
    labeler=x -> "", 
    kwargs...
)
    ax = Axis(
        fig[1, 1]; 
        xscale=log10, 
        xminorticksvisible=true, 
        xminorticks=IntervalsBetween(9),
        xticks=[0.1, 1.0, 10.0],
    )
    for (sim, c) in zip(sims, colors)
        viz(sim, config, fig, ax; color=c, labeler=labeler, kwargs...)
    end
    axislegend(ax)
    fig
end

function Utilities.viz(
    sim::WarrenTuningCurve, 
    config::Config, 
    fig::Figure, 
    ax::Axis; 
    color=:black,
    plot_cf=true,
    markersize=10.0,
    linewidth=2.0,
    correct_baseline=true,
    labeler=x->"",
    ylims=(-50.0, 0.0),
    xlims=(0.05, 40.0),
    kwargs...
)
    # Fetch data
    Δ = simulate(sim, config; kwargs...)

    # Correct for baseline, optionally
    if correct_baseline
        Δ .-= maximum(Δ)
    end

    # Plot
    scatter!(ax, Utilities.freq(sim) ./ 1e3, Δ; color=color, markersize=markersize)#, label=labeler(sim))
    lines!(ax, Utilities.freq(sim) ./ 1e3, Δ; color=color, linewidth=linewidth, label=labeler(sim))#, label=labeler(sim))
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)

    # Add lines and limits
    ax.xlabel = "Contralateral\nsuppressor freq (kHz)"
    ax.ylabel = "Rate suppression (%)"

    # If plot_cf, plot vertical line at cf
    if plot_cf
        vlines!(ax, getcf(sim)/1e3; color=color, linestyle=:dash)
    end

    fig
end