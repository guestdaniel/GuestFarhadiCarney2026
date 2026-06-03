export PrecursorRLF, PrecursorRLFSet, levels, construct_stimuli, precursors, probes, viz_breakout, viz_summary

"""
    PrecursorRLF

Simulation of RLF for probes following fixed precursor sound.

`simulate(sim)`
Runs the model on each acoustic waveform in the RLF separately. The returns are in the form 
of a single average rate for each stimulus within the specified analysis window provided by
idxswin(::PrecursorStimulus).

`construct_stimuli(sim)`
Construct and return each of the stimulus objects from each sound level for inspection.
"""
@with_kw struct PrecursorRLF{M, R, S} <: AbstractRLF where {M <: Model, R <: AbstractStimulus, S <: AbstractStimulus}
    model::M
    precursor::R
    probes::Vector{S}
    isi::Float64=0.3
    summaryfunc::Function=mean
end

function PrecursorRLF(
    precursor, 
    model; 
    levels_probe=0.0:5.0:80.0,
    dur_probe=0.01,
    dur_ramp_probe=0.005,
    isi=0.3,
    summaryfunc=mean,
)
    probes = map(levels_probe) do l
        PureTone(; 
            freq=first(model.cf), 
            dur=dur_probe, 
            dur_ramp=dur_ramp_probe, 
            fs=samprate(model), 
            level=l
        )
    end
    PrecursorRLF(; model=model, precursor=precursor, isi=isi, summaryfunc=summaryfunc, probes=probes)
end

Utilities.level(x::PrecursorRLF) = getfield.(x.probes, :level)

function Utilities.simulate(x::PrecursorRLF)
    @unpack_PrecursorRLF x
    # Simulate responses to each stimulus
    @info "Computing $(length(x.probes)) probe responses..."
    @showprogress pmap(construct_stimuli(x)) do stim
        # Compute response
        resp = compute(model, stim)

        # Compute mean rate in window of interest
        summaryfunc(resp[idxswin(stim)])
    end
end

function construct_stimuli(x::PrecursorRLF)
    @unpack_PrecursorRLF x
    map(probes) do probe
        PrecursorStimulus(; precursor=precursor, stimulus=probe, dur_isi=isi)
    end
end

"""
    PrecursorRLFSet

Simulation of series of related PrecursorRLFs .

Contains single field, `rlfs`, which is an M-tuple of N-tuples of RLFs. Here, M=number of 
different models tested, and N=number of different precursor stimuli tested.
"""
struct PrecursorRLFSet{M, N} <: Simulation where {M <: Int, N <: Int}
    rlfs::NTuple{M, NTuple{N, PrecursorRLF}}
end

function Base.hash(sim::PrecursorRLFSet)
    hashes = UInt64[]
    push!(hashes, Base.hash(typeof(sim)))
    map(sim.rlfs) do rlf_set
        map(rlf_set) do rlf
            push!(hashes, Base.hash(rlf))
        end
    end
    reduce(+, hashes)
end

function PrecursorRLFSet(; 
    cf=2.8e3,
    fs=100e3,
    params_model=[
        Dict(:moc_weight_wdr => 0.0, :moc_weight_ic => 0.0),
        Dict(:moc_weight_wdr => 5.0, :moc_weight_ic => 10.0),
    ],
    stage="hsr",
    species="human",
    level_pre=50.0,
    dur_pre=0.3,
    precursors=[
        PureTone(; freq=cf, dur=dur_pre, fs=fs, level=-Inf), 
        PureTone(; freq=cf, dur=dur_pre, fs=fs, level=level_pre),
    ],
    levels_probe=0.0:5.0:90.0,
    dur_probe=0.01,
    dur_ramp_probe=0.005,
    probes=[PureTone(; freq=cf, dur=dur_probe, dur_ramp=dur_ramp_probe, fs=fs, level=l) for l in levels_probe],
    isi=0.1,
    summaryfunc=mean,
)
    # Set up models
    models = [GFC2024(; cf=[cf], stage=stage, species=species, fs=fs, params...) for params in params_model]

    # Set up RLFs
    rlfs = map(models) do model 
        rlfs_this_model = map(precursors) do precursor
            PrecursorRLF(; model=model, precursor=precursor, probes=probes, isi=isi, summaryfunc=summaryfunc)
        end
        tuple(rlfs_this_model...)
    end
    PrecursorRLFSet(tuple(rlfs...))
end

function Utilities.simulate(sim::PrecursorRLFSet)
    map(x -> map(simulate, x), sim.rlfs)
end

Utilities.modelobj(sim::PrecursorRLFSet) = getfield.(first.(sim.rlfs), :model)
precursors(sim::PrecursorRLFSet) = getfield.(first(sim.rlfs), :precursor)
probes(sim::PrecursorRLFSet) = getfield(first(first(sim.rlfs)), :probes)
Utilities.level(sim::PrecursorRLFSet) = level(first(first(sim.rlfs)))
isi(sim::PrecursorRLFSet) = getfield(first(first(sim.rlfs)), :isi)

function Base.display(sim::PrecursorRLFSet)
    printstyled("$(typeof(sim))\n"; color=:blue, bold=true)
end

function Utilities.analyze(sim::PrecursorRLFSet; config=Config())
    # Summon data
    μs = @memo config simulate(sim)
    l = level(sim)

    # Calc ΔL for each pair of RLFs (first is assumed to reflect response to control 
    # stimulus, i.e., no precursor, subsequent RLFs are assumed to reflect response to
    # real precursor stimulus).
    map(μs) do μ_this_model
        baseline = fill(μ_this_model[1], length(μ_this_model[2:end]))
        map(zip(baseline, μ_this_model[2:end])) do (μ1, μ2)
            keys = ["ΔL", "ΔS", "CR", "L@CR", "idxs", "slopes", "intercepts"]
            values = calc_ΔL(l, μ1, μ2)
            Dict(keys .=> values)
        end
    end
end

function Utilities.viz(
    sim::PrecursorRLFSet; 
    config=Config(),
    colors=[:gray, colorschemes[:berlin]...],
    styles=[:dash, :solid, :dot]
)
    # Set up figure
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1:2, 1])

    # Summon data
    μs = @memo config simulate(sim)
    l = level(sim)

    # Create empty vectors to store plot elements
    markers = Scatter[]
    lines = Lines[]

    # Plot each
    map(enumerate(zip(μs, colors, modelobj(sim)))) do (idx, (μ_set, color, model))
        objs = map(zip(μ_set, styles, precursors(sim))) do (μ, style, precursor)
            # Plot
            _s = scatter!(ax, l, μ, color=color)
            _l = lines!(ax, l, μ, color=color, linestyle=style)

            return _s, _l
        end
        push!(markers, objs[1][1])
        if idx == 1
            [push!(lines, x[2]) for x in objs]
        end
    end

    # Add labels
    ax.xlabel = "Probe level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"
    ylims!(ax, 0.0, maximum(maximum.(maximum.(μs)))*1.2)
    ax.title = "ISI = $(round(isi(sim)*1000)) ms"

    # Create labels 
    lab_mod = map(model -> id(model; accesses=standardfields("MOC"), connector="  "), modelobj(sim))
    lab_pre = map(precursor -> id(precursor; accesses=[:level, :dur], connector = "  "), precursors(sim))

    Legend(fig[1, 2], markers, [lab_mod...], "Model")
    Legend(fig[2, 2], lines, [lab_pre...], "Precursor")

    fig
end

function viz_breakout(
    sim::PrecursorRLFSet; 
    config=Config(),
    n_precursor=length(precursors(sim)),
    colors=[:gray, get(colorschemes[:berlin], LinRange(0.0, 1.0, n_precursor))...],
    fig=Figure(size=(450*length(modelobj(sim)), 450)),
)
    # Determine number 
    n_model = length(modelobj(sim))
    n_precursor = length(precursors(sim))

    # Set up figure
    axs = [Axis(fig[2:3, i]) for i in 1:n_model]

    # Summon data
    μs = @memo config simulate(sim)
    l = level(sim)

    # Create empty vectors to store plot elements
    markers = Scatter[]
    lines = Lines[]

    # Loop over models and plot each's results
    map(enumerate(zip(μs, modelobj(sim), axs, sim.rlfs))) do (idx, (μ_set, model, ax, rlf_set))
        # Set up custom stimulus axes for this model 
        axs_stim = [Axis(fig[1, idx][i, 1]) for i in 1:n_precursor]
        hideydecorations!.(axs_stim)
        hidexdecorations!(axs_stim[end], ticks=false, ticklabels=false)
        hidexdecorations!.(axs_stim[1:(end-1)], ticks=false)

        # Loop over RLFs for each model and plot results
        objs = map(enumerate(zip(μ_set, colors, rlf_set))) do (idx_stim, (μ, color, rlf))
            # Plot stimulus example
            stims = construct_stimuli(rlf)
            w = synthesize(middle(stims))
            lines!(axs_stim[idx_stim], timevec(w, samprate(middle(stims))), w)

            # Plot rate-level function (μ vs l)
            _s = scatter!(ax, l, μ, color=color)
            _l = lines!(ax, l, μ, color=color)

            return _s, _l
        end
        push!(markers, objs[1][1])
        if idx == 1
            [push!(lines, x[2]) for x in objs]
        end
    end

    # Add labels
    [ax.xlabel = "Probe level (dB SPL)" for ax in axs]
    [ax.ylabel = "Firing rate (sp/s)" for ax in axs]
    ylims!.(axs, 0.0, maximum(maximum.(maximum.(μs)))*1.2)
    hideydecorations!.(axs[2:end], ticks=false, grid=false)

    # Create labels 
    lab_mod = map(model -> id(model; accesses=vcat(standardfields("MOC"), :stage), connector=" \n"), modelobj(sim))
    lab_pre = map(precursor -> id(precursor; accesses=[:level, :dur], connector = "  "), precursors(sim))

#    Legend(fig[2, length(modelobj(sim)) + 1], markers, [lab_mod...], "Model")
    [ax.title = x for (ax, x) in zip(axs, lab_mod)]
    Legend(fig[:, length(modelobj(sim)) + 1], lines, [lab_pre...], "Precursor")

    fig
end

function viz_summary(
    x::AbstractVector,
    sim::PrecursorRLFSet; 
    xlabel="Spectrum level (dB SPL)",
    config=Config(),
    colors=[:gray, colorschemes[:Dark2_8]...],
    fig=Figure(size=(450, 450)),
    ax_ΔL=Axis(fig[1, 1]),
)
    # Fetch results
    results = analyze(sim; config=config)
    ΔL = map(x -> getindex.(x, "ΔL"), results)

    # Plot ΔL results
    lns = map(zip(ΔL, colors)) do (_ΔL, color)
        linescatter!(ax_ΔL, x, _ΔL; color=color)
    end
    xlims!(ax_ΔL, extrema(x)...)
    ylims!(ax_ΔL, 0.0, 25.0)
    ax_ΔL.xlabel = xlabel
    ax_ΔL.ylabel = "ΔL (dB)"
    cf = first(getfield(first(modelobj(sim)), :cf))
    species = getfield(first(modelobj(sim)), :species)
    guardrail = gr_ΔL_trendline(f_to_cd(cf, species), gr_optimize_ncd_gaussian()...)
    hlines!(ax_ΔL, guardrail; color=:gray, linestyle=:dash)

    # Prepare legend
    labels = id.(modelobj(sim); accesses=standardfields("MOC"), connector=" ")
#    Legend(fig[1, 2], lns, labels)

    # Return and display figure
    display(fig)
    fig, ax_ΔL
end
