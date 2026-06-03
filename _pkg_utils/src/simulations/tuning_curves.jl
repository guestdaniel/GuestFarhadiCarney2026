export AbstractTuningCurve
export RLFTC
export IsoLevelTC
export threshold_curve

abstract type AbstractTuningCurve <: Simulation end

freq(sim::AbstractTuningCurve) = sim.freqs

"""
    RLFTC <: AbstractTuningCurve

Tuning curve derived from rate-level functions

Tuning curve estimated from rate-level functions (RLFs) at different frequencies. 
The tradition tuning curve (threshold vs frequency) is computed by `threshold_curve`, which
in turn simply calls `threshold` on each RLF object.

# Fields
- `model::M`: Model object 
- `rlfs::Vector{S}`: Vector of RLF objects
- `summaryfunc::Function`: Function to summarize responses
"""
struct RLFTC{M, S} <: AbstractTuningCurve where {M <: Model, S <: AbstractRLF}
    model::M
    rlfs::Vector{S}
    summaryfunc::Function
end

function RLFTC(
    model::Model; 
    freq_low=-1.0,
    freq_high=1.0,
    n_freq=21,
    freqs=OctRange(getcf(model), freq_low, freq_high, n_freq), 
    level_low=0.0,
    level_step=2.5,
    level_high=70.0,
    levels=level_low:level_step:level_high,
    kwargs...
)
    rlfs = map(freqs) do freq
        RLF(model; freq=freq, levels=levels, kwargs...)
    end
    RLFTC(model, rlfs, mean)
end

level(sim::RLFTC) = level(sim.rlfs[1])

freq(sim::RLFTC) = freq.(sim.rlfs)

function simulate(sim::RLFTC)
    map(simulate, sim.rlfs)
end

function threshold_curve(sim::RLFTC, cutoff=0.1; config=Config())
    rlfs = @memo config simulate(sim)
    l = level(sim)
    map(x -> threshold(l, x, cutoff), rlfs)
end

function threshold(sim::RLFTC; config=Config())
    minimum(threshold_curve(sim; config=config))
end

function analyze(sim::RLFTC; config=Config())
    analyze_tc(freq(sim), threshold(sim; config=config))
end

function viz(
    sim::RLFTC; 
    fig=Figure(), 
    ax=Axis(fig[1, 1]; xscale=log10),
    cutoff=0.1,
    config=Config(),
)
    if typeof(cutoff) <: AbstractFloat
        lines!(ax, freq(sim) ./ 1e3, threshold_curve(sim, cutoff; config=config))
    else
        map(cutoff) do _c
            lines!(ax, freq(sim) ./ 1e3, threshold_curve(sim, _c; config=config))
        end
    end
    xlims!(ax, 0.1, 32.0)
    ylims!(ax, -10.0, 90.0)
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xminorticksvisible = true
    ax.xminorticks = IntervalsBetween(9)
    display(fig)
    return fig, ax
end

"""
    analyze_tc(freqs, levels)

Analyzes an iso-response tuning curve, returning CF, threshold, Q10 and Q20
"""
function analyze_tc(freqs, levels)
    # Extract CF and threshold
    idx_cf = argmin(levels)
    cf = freqs[idx_cf]
    θ = levels[idx_cf]

    # Extract Q values
    Qs = map([10.0]) do shift
        # Create temporary variable levels .- θ .- shift, which is zero at our points of
        # interest (where tuning curve is 10 dB above threshold)
        temp = levels .- θ .- shift

        # Find indices closest to zero-crossings (leading index)
        idxs = findall(diff(sign.(temp)) .!= 0.0)
        if length(idxs) < 2
            return NaN
        end

        # If freqs[idxs] don't straddle CF, return early with NaN
        if (freqs[idxs[1]] > cf) | (freqs[idxs[2]] < cf)
            return NaN
        end

        # Determine frequency where curve crosses zero by linear interpolation on 
        # log-frequency axis
        f_intersect = map(idxs) do idx
            β₁ = (temp[idx+1] - temp[idx]) / log2(freqs[idx+1] / freqs[idx])
            β₀ = temp[idx] - β₁ * log2(freqs[idx]/cf)
            cf * 2.0 ^ (-β₀/β₁)
        end

        # Calculate bandwidth as difference of two crossing points
        bw = f_intersect[2] - f_intersect[1]

        # Return Q
        return cf/bw
    end

    # Return
    return cf, θ, Qs...
end

function analyze_tc(freqs, levels...)
    map(l -> analyze_tc(freqs, l), levels)
end

"""
    IsoLevelTC <: AbstractTuningCurve

Iso-level tuning curve 

Tuning curve consisting of response at different frequencies at a fixed level.

# Fields
- `model::M`: Model object 
- `probes::Vector{S}`: Vector of probe tones
- `summaryfunc::Function`: Function to summarize responses
"""
struct IsoLevelTC{M, S} <: AbstractTuningCurve where {M <: Model, S <: AbstractStimulus}
    model::M
    probes::Vector{S}
    summaryfunc::Function
end

function simulate(sim::IsoLevelTC)
    map(sim.probes) do probe
        r = compute(sim.model, probe)
        sim.summaryfunc(r[idxswin(probe)])
    end
end

function freq(sim::IsoLevelTC)
    freq.(sim.probes)
end