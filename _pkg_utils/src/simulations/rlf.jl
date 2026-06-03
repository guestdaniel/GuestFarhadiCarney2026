export AbstractRLF, RLF, plot_rlf, midpoint, threshold, dynamicrange, threshold_shift, midpoint_shift, raterange

abstract type AbstractRLF <: Simulation end

"""
    RLF

Simulate responses to probe stimuli at different levels for given model
"""
struct RLF{M, S} <: AbstractRLF where {M <: Model, S <: AbstractStimulus}
    model::M
    probes::Vector{S}
    summaryfunc::Function
end

function RLF(
    model::M; 
    fs=100e3, 
    level_min=-5.0,
    level_max=80.0,
    level_step=5.0,
    levels=level_min:level_step:level_max,
    dur=0.1,
    summaryfunc=mean, 
    kwargs...
) where {M <: Model}
    # Create Stimulus objects
    probes = map(l -> PureTone(; freq=isa(getcf(model), Vector) ? middle(getcf(model)) : getcf(model), level=l, fs=fs, dur=dur, kwargs...), levels)

    # Create RLF
    RLF(model, probes, summaryfunc)
end

level(sim::AbstractRLF) = level.(sim.probes)
freq(sim::AbstractRLF) = freq(sim.probes[1])
stimobj(sim::AbstractRLF) = sim.probes

function simulate(sim::AbstractRLF) 
    pmap(sim.probes) do probe
        # Simulate response
        r = compute(sim.model, probe)

        # Extract time window
        sim.summaryfunc(r[idxswin(probe)])
    end
end

"""
    smooth(sim::AbstractRLF, x, y)

Smooths output of sim(rlf) via fitting logistic function (needs hard saturation)
"""
function smooth(::AbstractRLF, x, y)
    # Guestimate reasonable parameter values to start with
    r_min = minimum(y)
    r_max = maximum(y)
    rr = raterange(y)
    θ₃ = r_min
    θ₄ = r_max
    θ₁ = x[argmin(abs.(y .- (r_min .+ rr/2)))]

    # Pass to smooth and return x̂ and ŷ
    smooth(x, y, logistic_good; x₀=[θ₁, 0.25, θ₃, θ₄])
end

"""
    raterange(sim::AbstractRLF[; config=Config()])
"""
function raterange(sim::AbstractRLF; config=Config(), interp=false)
    # Fetch RLF data
    μ = @memo config simulate(sim)

    # Optionally interpolate and then extract rate range
    if interp
        _, μ = smooth(sim, level(sim), μ)
    end

    raterange(μ)
end

function raterange(μ::AbstractVector)
    maximum(μ) - minimum(μ)  # difference between highest and lowest rate
end

"""
    midpoint(sim::AbstractRLF[; config=Config()])

Compute midpoint of deterministic RLF as level closest to achieving 50% of total rate range above spont
"""
function midpoint(sim::AbstractRLF; config=Config(), interp=false)
    # Fetch RLF data
    μ = @memo config simulate(sim)

    # Optionally interpolate and then extract midpoint
    l = level(sim)
    if interp
        l, μ = smooth(sim, l, μ)
    end

    nthpoint(l, μ, 0.5)
end

function nthpoint(l::AbstractVector, μ::AbstractVector, p::Float64)
    ΔR = raterange(μ)
    l[argmin(abs.(μ .- (minimum(μ) + ΔR*p)))]
end

"""
    threshold(sim::AbstractRLF[; config=Config()])

Compute threshold for deterministic RLF as level achieving 10% of rate range above spont
"""
function threshold(sim::AbstractRLF, cutoff=0.1; config=Config(), interp=false)
    # Fetch RLF data
    μ = @memo config simulate(sim) 

    # Optionally interpolate and then extract threshold
    l = level(sim)
    if interp
        l, μ = smooth(sim, l, μ)
    end

    threshold(l, μ, cutoff)
end

function threshold(l::AbstractVector, μ::AbstractVector, cutoff=0.1)
    # Determine rate range of RLF
    ΔR = raterange(μ)

    # If ΔR is less than cutoff*spont return NaN
    if ΔR < cutoff*minimum(μ)
        return NaN
    end

    # Find levels that are greater than spont+cutoff*RR 
    idxs = findall(μ .> (minimum(μ) + ΔR*cutoff))

    # If we didn't find any, return NaN, otherwise return level at first hit index
    if isempty(idxs)
        return NaN
    else
        return l[first(idxs)]
    end
end

"""
    dynamicrange(sim::AbstractRLF; config=Config())

Compute dynamic range for deterministic RLF as span from spont+10% of rate range to spont+90% of level range; return DR value in dB and bottom and top of range in dB.
"""
function dynamicrange(sim::AbstractRLF; config=Config(), interp=false)
    μ = @memo config simulate(sim) 
    l = level(sim)
    if interp
        l, μ = smooth(sim, l, μ)
    end
    dynamicrange(l, μ)
end

function dynamicrange(l::AbstractVector, μ::AbstractVector)
    # Compute rate range 
    ΔR = raterange(μ)

    # Determine levels closest to yielding response at min+0.1ΔR to min+0.9ΔR (i.e., max-0.1ΔR)
    idxs_bottom = findall(μ .> (minimum(μ) + 0.1ΔR))
    idxs_top = findall(μ .< (maximum(μ) - 0.1ΔR))
    if ~(isempty(idxs_bottom) & isempty(idxs_top))
        idx_bottom = first(idxs_bottom)
        idx_top = last(idxs_top)
        level_bottom = l[idx_bottom]
        level_top = l[idx_top]
        dr = level_top - level_bottom
    else
        dr, level_bottom, level_top, idx_bottom, idx_top = NaN, NaN, NaN, NaN, NaN
    end
    
    # Return dynamic range, levels, and indices
    return dr, (level_bottom, level_top), (idx_bottom, idx_top)
end

"""
    threshold_shift(rlf1::AbstractRLF, rlf2::AbstractRLF[; config=Config()])

Compute threshold shift for two RLFs as simple difference of their thresholds
"""
function threshold_shift(rlf1::AbstractRLF, rlf2::AbstractRLF; config=Config(), interp=false)
    # Fetch or simultate RLFs
    μ1 = @memo config simulate(rlf1)
    μ2 = @memo config simulate(rlf2)

    # Determine levels used
    l = level(rlf1)
    @assert all(l .== level(rlf2))  # check to make sure RLFs used same levels

    # Optionally interpolate the RLFs and then threshold shift
    if interp
        l̂, μ̂2 = smooth(rlf1, l, μ2)
        l̂, μ̂1 = smooth(rlf2, l, μ1)
        threshold(l̂, μ̂2) - threshold(l̂, μ̂1)
    else
        threshold(l, μ2) - threshold(l, μ1)
    end
end

"""
    midpoint_shift(rlf1::AbstractRLF, rlf2::AbstractRLF[; config=Config()])

Compute midpoint shift for two RLFs as simple difference of their midpoints
"""
function midpoint_shift(rlf1::AbstractRLF, rlf2::AbstractRLF; config=Config(), interp=false)
    # Fetch or simultate RLFs
    μ1 = @memo config simulate(rlf1)
    μ2 = @memo config simulate(rlf2)

    # Determine levels used
    l = level(rlf1)
    @assert all(l .== level(rlf2))  # check to make sure RLFs used same levels

    # Optionally interpolate the RLFs and then estimate thresholds
    if interp
        l̂, μ̂2 = smooth(rlf2, l, μ2)
        l̂, μ̂1 = smooth(rlf1, l, μ1)
        nthpoint(l̂, μ̂2, 0.5) - nthpoint(l̂, μ̂1, 0.5)
    else
        nthpoint(l, μ2, 0.5) - nthpoint(l, μ1, 0.5)
    end
end

# Define some visualization methods for rate-level functions
function viz(sim::AbstractRLF; config=Config(), interp=false, kwargs...)
    μ = @memo config simulate(sim)
    if interp
        fig, ax = plot_rlf(smooth(sim, level(sim), μ)...; kwargs...)
        scatter!(ax, level(sim), μ; markersize=10.0)
    else
        fig, ax = plot_rlf(level(sim), μ; kwargs...)
    end
    display(fig)
    fig, ax
end

function viz(sims::Vector{<:AbstractRLF}; config=Config(), kwargs...)
    μ = map(sims) do sim 
        @memo config simulate(sim)
    end
    plot_rlf(level(sims[1]), μ; kwargs...)
end

function plot_rlf(
    l::AbstractVector, 
    μ::AbstractVector; 
    fig=Figure(), 
    yscale=identity,
    ax=Axis(fig[1, 1]; yscale=yscale), 
    color=:black,
    label="",
    linestyle=:solid,
    linewidth=1.0,
    markersize=6.0,
    plot_ΔL=true,
    ΔL=5.0,
    ref_levels=[20.0, 45.0],
    ref_level_colors=colorschemes[:Dark2_8],
    ylims=nothing,
    ylim_mode="auto",
    xlims=nothing,
    xlabel="Probe level (dB SPL)",
    ylabel="Firing rate (sp/s)",
)
    # Plot curve
    if ax.yscale[] == log10
        μ = max.(μ, 0.0) .+ 1e-12  # avoid log10(0) issues
    end

    # If plotting ΔL, draw lines connecting reference levels to incremented levels
    for (ref_level, c) in zip(ref_levels, ref_level_colors)
        if plot_ΔL && any(l .== ref_level) && any(l .== ref_level + ΔL)
            idx_ref = findfirst(==(ref_level), l)
            idx_inc = findfirst(==(ref_level + ΔL), l)
            lines!(ax, [ref_level, ref_level, maximum(l)], [0.0, μ[idx_ref], μ[idx_ref]]; color=(c, 0.5), linestyle=:dash, linewidth=linewidth/2)
            lines!(ax, [ref_level + ΔL, ref_level + ΔL, maximum(l)], [0.0, μ[idx_inc], μ[idx_inc]]; color=(c, 0.5), linestyle=:dash, linewidth=linewidth/2)
        end
    end

    # Plot data
    lines!(ax, l, μ, color=color, linestyle=linestyle, linewidth=linewidth)

    # Add labels
    ax.xlabel = xlabel
    ax.ylabel = ylabel

    # Set ylabel to be in dB if using log scale
    if yscale != identity
        ticks = maximum(μ) .* 10 .^ (-100.0:10.0:0.0 ./ 20) 
        ax.yticks = ticks
    end
    ax.xticks = floor(minimum(l); digits=-1):20.0:ceil(maximum(l); digits=-1)
    ax.xminorticksvisible = true
    ax.xminorticks = floor(minimum(l); digits=-1):5.0:ceil(maximum(l); digits=-1)

    # Set limits
    if ylims !== nothing
        ylims!(ax, ylims...)
    else
        ylims!(ax, ax.yscale[] == log10 ? minimum(μ) / 1.2 : 0.0, ylim_mode == "auto" ? maximum(μ)*1.2 : 1000.0)
    end
    if xlims !== nothing
        xlims!(ax, xlims...)
    end
    display(fig)
    fig, ax
end

function plot_rlf(
    l::Vector{Float64}, 
    μ::Vector{Vector{Float64}}; 
    fig = Figure(),
    yscale=identity,
    ax = Axis(fig[1, 1]; yscale=yscale),
    colors=colorschemes[:Dark2_8],
    labels=fill("", length(colors)),
    show_legend=true,
    position=:rb,
    kwargs...
)
    map(zip(μ, colors, labels)) do (_μ, c, lab)
        plot_rlf(l, _μ; fig=fig, ax=ax, color=c, label=lab, kwargs...)
    end
    if yscale == identity
        ylims!(ax, 0.0, maximum(maximum.(μ)) * 1.2)
    end
    # if show_legend axislegend(; position=position) end
    display(fig)
    fig, ax
end