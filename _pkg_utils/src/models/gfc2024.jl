export GFC2024, mocparams, moc_nonlin, plot_moc_nonlin

"""
    GFC2024 <: Model

Auditory efferent model of Guest, Farhadi, and Carney (2024) 
"""
@with_kw struct GFC2024 <: Model
    # General parameters
    fs::Float64 = 100e3
    cf::Vector{Float64} = [1000.0]
    coi::Vector{Int64} = [Int(ceil(length(cf) / 2))]
    n_chan::Int64 = length(cf)
    cf_low::Float64 = minimum(cf)
    cf_high::Float64 = maximum(cf)
    audiogram::Audiogram = Audiogram()
    stage::String = "hsr"

    # IHC parameters
    cohc::Vector{Float64} = audiogram.desc == "NH" ? ones(size(cf)) : fit_audiogram(audiogram, cf)[1]
    cihc::Vector{Float64} = audiogram.desc == "NH" ? ones(size(cf)) : fit_audiogram(audiogram, cf)[2]
    species::String = "human"

    # ANF parameters
    powerlaw_mode::Int64 = 2
    fractional::Bool = false

    # IC parameters
    cn_tau_e::Float64 = 0.5e-3
    cn_tau_i::Float64 = 2.0e-3
    cn_delay::Float64 = 1.0e-3
    cn_amp::Float64 = 1.5
    cn_inh::Float64 = 0.6
    ic_tau_e::Float64 = 1.0 / (10.0 * 64.0)
    ic_tau_i::Float64 = ic_tau_e * 1.5
    ic_delay::Float64 = ic_tau_e * 2.0
    ic_amp::Float64 = 1.0
    ic_inh::Float64 = 0.9

    # MOC parameters
    moc_cutoff::Float64 = 0.64
    moc_beta::Vector{Float64} = fill(0.015, length(cf))
    moc_offset::Vector{Float64} = fill(0.0, length(cf))
    moc_minval::Float64 = 0.0
    moc_maxval::Float64 = 1.0
    moc_weight::Vector{Float64} = fill(0.0, length(cf))
    moc_width::Float64 = 1.0
    moc_delay::Float64 = 0.025

    # Other params
    dur_pad_left::Float64 = 0.02  # also used as "dur_settle"
    clip_left::Bool = dur_pad_left == 0.0 ? false : true
    dur_pad_right::Float64 = 0.02
    clip_right::Bool = false
end

# Some constructors for convenience
GFC2024(cf::Float64; kwargs...) = GFC2024(; cf=[cf], kwargs...)
GFC2024(cf::Vector{Float64}; kwargs...) = GFC2024(; cf=cf, kwargs...)

# Some other functions extending
islogout(m::GFC2024) = any(m.stage .== ["c1", "c2", "ihc"])

# Function to convert from the strings used to indicate stages in the dict model functions to fieldnames used in the memory structure
function stage_to_fieldname(stage::String)
    if stage == "hsr"
        return :hsrout
    elseif stage == "lsr"
        return :lsrout
    elseif stage == "ihc"
        return :ihcout
    elseif stage == "c1"
        return :c1out
    elseif stage == "c2"
        return :c2out
    elseif stage == "control"
        return :controlout
    elseif stage == "cn"
        return :cnout
    elseif stage == "ic"
        return :icout
    elseif stage == "expon_hsr"
        return :expout_hsr
    elseif stage == "syn_hsr"
        return :synout_hsr
    elseif stage == "expon_lsr"
        return :expout_lsr
    elseif stage == "syn_lsr"
        return :synout_lsr
    else
        return Symbol(stage)
    end
end

function stage_to_index(stage::String)
    if stage == "control"
        return 1
    elseif stage == "c1"
        return 2
    elseif stage == "c2"
        return 3
    elseif stage == "ihc"
        return 4
    elseif stage == "expon_hsr"
        return 5
    elseif stage == "sout1_hsr"
        return 6
    elseif stage == "sout2_hsr"
        return 7
    elseif stage == "syn_hsr"
        return 8
    elseif stage == "expon_lsr"
        return 9
    elseif stage == "sout1_lsr"
        return 10
    elseif stage == "sout2_lsr"
        return 11
    elseif stage == "syn_lsr"
        return 12
    elseif stage == "hsr"
        return 13
    elseif stage == "lsr"
        return 14
    elseif stage == "cn"
        return 15
    elseif stage == "ic"
        return 16
    elseif stage == "mocwdr"
        return 17
    elseif stage == "mocic"
        return 18
    elseif stage == "gain"
        return 19
    elseif stage == "gainpostmix"
        return 20
    else
        error("Unknown stage: $stage")
    end
end

# compute(model, stimulus) maps from stimulus to response
function _compute(m::GFC2024, x::Vector{Float64})
    sim_gfc2023_dict(
        x,
        m.cf;
        # General parameters
        fs=m.fs,
        # IHC parameters
        cohc=m.cohc,
        cihc=m.cihc,
        species=m.species,
        # ANF parameters
        fractional=m.fractional,
        powerlaw_mode=m.powerlaw_mode,
        # IC parameters
        cn_tau_e=m.cn_tau_e,
        cn_tau_i=m.cn_tau_i,
        cn_delay=m.cn_delay,
        cn_amp=m.cn_amp,
        cn_inh=m.cn_inh,
        ic_tau_e=m.ic_tau_e,
        ic_tau_i=m.ic_tau_i,
        ic_delay=m.ic_delay,
        ic_amp=m.ic_amp,
        ic_inh=m.ic_inh,
        # MOC parameters
        moc_cutoff=m.moc_cutoff,
        moc_beta=m.moc_beta,
        moc_offset=m.moc_offset,
        moc_minval=m.moc_minval,
        moc_maxval=m.moc_maxval,
        moc_weight=m.moc_weight,
        moc_width=m.moc_width,
        moc_delay=m.moc_delay,
        # Misc parameters
        dur_pad_left=m.dur_pad_left,
        clip_left=m.clip_left,
        dur_pad_right=m.dur_pad_right,
        clip_right=m.clip_right,
    )
end

# Compute methods mapping from stimulus to response
_compute(m::GFC2024, x::AbstractStimulus) = _compute(m, synthesize(x))
compute(m::GFC2024, x::AbstractBinauralStimulus) = compute(m, ipsi(x))
_compute(m::GFC2024, x::AbstractBinauralStimulus) = _compute(m, synthesize(ipsi(x)))

function compute(m::GFC2024, x::Vector{Float64})
    extract(m, _compute(m, x)[m.stage])
end

# Compute method with reused memory
function compute(m::GFC2024, x::AbstractStimulus, N::Int64)
    # Synthesize x N times
    xs = [synthesize(x) for i in 1:N]

    # Create memory
    mem = GFC2023_Mem(xs[1], m.cf; fs=samprate(m), dur_pad_left=m.dur_pad_left, dur_pad_right=m.dur_pad_right)

    # Pass to specific _compute method
    _compute(m, xs, mem)
end

function _compute(m::GFC2024, xs::Vector{Vector{Float64}}, mem::GFC2023_Mem)
    # Accumulate and extract
    map(xs) do x
        # Compute model response in-place and extract
        sim_gfc2023!(
            mem,
            x,
            m.cf;
            # IHC parameters
            cohc=m.cohc,
            cihc=m.cihc,
            species=m.species,
            # ANF parameters
            fractional=m.fractional,
            powerlaw_mode=m.powerlaw_mode,
            # IC parameters
            cn_tau_e=m.cn_tau_e,
            cn_tau_i=m.cn_tau_i,
            cn_delay=m.cn_delay,
            cn_amp=m.cn_amp,
            cn_inh=m.cn_inh,
            ic_tau_e=m.ic_tau_e,
            ic_tau_i=m.ic_tau_i,
            ic_delay=m.ic_delay,
            ic_amp=m.ic_amp,
            ic_inh=m.ic_inh,
            # MOC parameters
            moc_cutoff=m.moc_cutoff,
            moc_beta=m.moc_beta,
            moc_offset=m.moc_offset,
            moc_minval=m.moc_minval,
            moc_maxval=m.moc_maxval,
            moc_weight=m.moc_weight,
            moc_width=m.moc_width,
            moc_delay=m.moc_delay,
            # Misc parameters
            dur_pad_left=m.dur_pad_left,
            clip_left=m.clip_left,
            dur_pad_right=m.dur_pad_right,
            clip_right=m.clip_right,
        )[stage_to_index(m.stage)]
    end
end


"""
    moc_nonlin(x[; w=1.0, β=0.015, θ=40.0, minval=0.1, maxval=1.0])

Compute output of pointwise MOC input-output nonlinearity

Computes output of rational input-output nonlinearity between LSR rates and output
gain factor, given weight on input `w`, slope factor `β`, threshold w.r.t. to weighted rate 
`θ`.
"""
function moc_nonlin(x; w=1.0, β=0.015, θ=40.0, minval=0.1, maxval=1.0)
    if (w * x) < θ
        return maxval
    else
        return ((maxval - minval) * 1.0 / (1.0 + (β * ((w * x) - θ))^2)) + minval
    end
end

function mocparams(model::GFC2024)
    @unpack_GFC2024 model
    # Print title
    printstyled("Efferent AN model (2024)\n"; color=:blue, bold=true)
    printstyled("$n_chan CFs [$(round(cf_low/1000; digits=2)) to $(round(cf_high/1000; digits=2)) kHz] at $(fs/1000) kHz, $stage\n")

    printstyled("MOC stage: ", bold=true, italic=true)
    print("delay = $(round(moc_delay*1000; digits=1)) ms\n")

    printstyled("MOC-WDR stage: ", bold=true, italic=true)
    params = round.([moc_weight, moc_beta, moc_offset, moc_maxval]; digits=3)
    strs = map(zip(["weight", "β", "TH", "max"], params)) do (name, param)
        "$name = $param"
    end
    print("$(join(strs, ", "))\n")
end

function plot_moc_nonlin(
    model::GFC2024;
    fig=Figure(; size=(450, 350)),
    ax=Axis(fig[1, 1]; yscale=log10),
    color=:black,
    label="",
)
    @unpack_GFC2024 model
    rates = 1:1:500
    lines!(ax, rates, moc_nonlin.(rates; w=moc_weight[coi][1], β=moc_beta, θ=moc_offset, minval=moc_minval, maxval=moc_maxval), color=color, label=label)
    ax.yticks = 0.0:0.2:1.0
    ax.xticks = 0.0:20.0:100.0
    xlims!(ax, 0.0, 75.0)
    ylims!(ax, 0.15, 1.1)
    ax.xlabel = "MOC rate (sp/s)"
    ax.ylabel = "Gain factor"
    fig
end

function viz(
    model::GFC2024,
    x::Vector{Float64};
    stages=["ihc", "c1", "c2", "control", "hsr", "lsr", "ic", "mocwdr", "mocic", "gain", "gainpostmix"],
)
    # Compute response
    resp = _compute(model, x)

    # Create figure and plot each element
    fig = Figure(; size=(600, max(500, 100 * length(stages))))
    axs = [Axis(fig[i, 1]) for i in 1:(length(stages)+1)]
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    lines!(axs[1], timevec(x, samprate(model)), x)
    for (idx, key) in enumerate(stages)
        r = resp[key][model.coi][1]
        lines!(axs[idx+1], timevec(r, samprate(model)), r)
    end
    xlims!.(axs, (extrema(timevec(x, samprate(model))) .+ (0.0, 0.025))...)
    for (idx, name) in enumerate(vcat("Stimulus", stages))
        Label(fig[idx, 2], name; tellheight=false)
    end
    fig
end

function viz(
    model::GFC2024,
    x::Vector{Vector{Float64}};
    stages=["ihc", "c1", "c2", "control", "sout1", "sout2", "expon", "syn", "hsr", "lsr", "ic", "mocwdr", "mocic", "gain", "gainpostmix"],
    colors=colorschemes[:Dark2_8],
)
    # Compute response
    resps = map(x) do _x
        _compute(model, _x)
    end

    # Create figure and plot each element
    fig = Figure(; size=(600, max(500, 100 * length(stages))))
    axs = [Axis(fig[i, 1]) for i in 1:(length(stages)+1)]
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    map(zip(x, colors)) do (_x, c)
        lines!(axs[1], timevec(_x, samprate(model)), _x; color=c)
    end
    for (idx, key) in enumerate(stages)
        map(zip(resps, colors)) do (resp, c)
            r = resp[key][model.coi][1]
            lines!(axs[idx+1], timevec(r, samprate(model)), r; color=c)
        end
    end
    xlims!.(axs, (extrema(vcat(timevec.(x, samprate(model))...)) .+ (0.0, 0.025))...)
    for (idx, name) in enumerate(vcat("Stimulus", stages))
        Label(fig[idx, 2], name; tellheight=false)
    end
    fig
end

viz(m::GFC2024, s::AbstractStimulus; kwargs...) = viz(m, synthesize(s); kwargs...)
viz(m::GFC2024, s::AbstractBinauralStimulus; kwargs...) = viz(m, synthesize(ipsi(s)); kwargs...)
viz(m::GFC2024, s::Vector{<:AbstractStimulus}; kwargs...) = viz(m, synthesize.(s); kwargs...)

function neurogram(
    model::Model,
    resp::Dict;
    stages=["hsr", "ic", "gain"],
    scale=1.25,
    size=(Int(round(750 * scale)), Int(round(300 * scale)) * length(stages)),
    config=Config(),
    colorrange_gain=(1e-1, 1.0),
    kwargs...
)
    # Set up figure
    fig = Figure(; size=size)

    # Loop through requested stages and build plot
    for (idx, stage) in enumerate(stages)
        # Create axis
        ax = Axis(fig[idx, 1]; yscale=log10)

        # Customize according to stage and plot neurogram
        if stage == "gain" || stage == "gainpostmix"
            # Plot neurogram with batlow and standard colorbar
            hm = neurogram!(ax, model.cf, resp[stage]; colormap=Reverse(colorschemes[:batlow]), colorrange=colorrange_gain, colorscale=log10, kwargs...)
            Colorbar(fig[idx, 2], hm; label=uppercase(stage))

            # Add magic RLF legend
            ax = Axis(fig[idx, 3][2, 1])
            cohc = LogRange(0.1, 1.0, 9)
            colors = get(reverse(colorschemes[:batlow]), cohc)
            map(zip(cohc, colors)) do (cohc, color)
                sim = RLF(GFC2024(2e3; moc_weight=[0.0], stage="hsr", cohc=[cohc]))
                μ = @memo config simulate(sim)
                lines!(ax, level(sim), μ; color=color, linewidth=2.0)
            end
            ylims!(ax, 0.0, 350.0)
            ax.ylabel = "Firing rate (sp/s)"
            ax.xlabel = "Probe level (dB SPL)"
            colsize!(fig.layout, 3, Relative(0.2))
        else
            # Plot neurogram with viridis and use standard colorbar
            hm = neurogram!(ax, model.cf, resp[stage]; colormap=colorschemes[:viridis], colorrange=standardylims(stage), kwargs...)
            Colorbar(fig[idx, 2], hm; label=uppercase(stage))
        end
    end
    fig
end

# Stimulus method, uses @memo internally to cache output of _compute, allowing reuse 
function neurogram(m::Model, s::AbstractStimulus; kwargs...)
    resp = _compute(m, s)
    neurogram(m, resp; kwargs...)
end

# Final method, fallback for arbitrary time-pressure waveform
function neurogram(model::Model, x::Vector{Float64}; kwargs...)
    resp = _compute(model, x)
    neurogram(model, resp; kwargs...)
end