export SLF, magnitude, RLFPair

# Below, we extend definitions of RLF constructors provided in Utilities.jl to use
# the ProbePureToneElicitor class to allow simulation of CAS RLFs
function Utilities.RLF(
    model::Model, 
    type::Type{ProbePureToneElicitor3}; 
    fs=100e3, 
    level_step=2.0,
    levels=0.0:level_step:60.0, 
    level_elicitor=70.0,
    freq_probe=getcf(model),
    freq_elicitor=getcf(model),
    summaryfunc=mean,
    kwargs...
)
    # Create Stimulus objects
    probes = map(l -> ProbePureToneElicitor3([freq_probe, l], [freq_elicitor, level_elicitor]; fs=fs, kwargs...), levels)

    # Create RLF
    Utilities.RLF(model, probes, summaryfunc)
end

function Utilities.RLF(
    model1::Model, 
    model2::Model,
    type::Type{ProbePureToneElicitor3}; 
    kwargs...
)
    # Assert that the two models have the same CFs and if so create RLFs
    @assert getcf(model1) == getcf(model2) 
    RLF(model1, ProbePureToneElicitor3; kwargs...), RLF(model2, ProbePureToneElicitor3; kwargs...)
end

function RLFPair(
    model::Model, 
    levels_elicitor::Vector{Float64}; 
    kwargs...
)
    # Assert that the two models have the same CFs and if so create RLFs
    RLF(model, ProbePureToneElicitor3; level_elicitor=levels_elicitor[1], kwargs...), 
    RLF(model, ProbePureToneElicitor3; level_elicitor=levels_elicitor[2], kwargs...)
end


# Here, we provide a constructor for an IsoLevelTC defined for the ProbePureToneElicitor 
# stimulus.
function Utilities.IsoLevelTC(
    model::Model, 
    type::Type{ProbePureToneElicitor3}; 
    fs=100e3, 
    freq_probe=getcf(model),
    level_probe=20.0,
    freqs_elicitor=OctRange(getcf(model), -3.0, 3.0, 21),
    level_elicitor=70.0,
    summaryfunc=mean,
    kwargs...
)
    # Create Stimulus objects
    probes = map(f -> ProbePureToneElicitor3([freq_probe, level_probe], [f, level_elicitor]; fs=fs, kwargs...), freqs_elicitor)

    # Create RLF
    Utilities.IsoLevelTC(model, probes, summaryfunc)
end

Utilities.freq(sim::IsoLevelTC{M, ProbePureToneElicitor3}) where {M <: Model} = map(x -> Utilities.freq(contra(x)), sim.probes)


# Here, we provide a definition for a suppression-level function (SLF), which is an 
# abstraction useful for determining magnitudes and thresholds for efferent suppression.
struct SLF{M, S} <: AbstractRLF where {M <: Model, S <: AbstractStimulus}
    model::M
    probes::Vector{S}
    summaryfunc::Function
end

function SLF(
    model::Model; 
    fs=100e3, 
    level_probe=20.0,
    level_elicitor_min=0.0,
    level_elicitor_max=80.0,
    level_elicitor_step=5.0,
    levels_elicitor=level_elicitor_min:level_elicitor_step:level_elicitor_max,
    kwargs...
)
    # Create probes + elicitors (ProbePureToneElicitors) with varying elicitor levels
    cf = getcf(model)
    probes = map(l -> ProbePureToneElicitor3([cf, level_probe], [cf, l]; fs=fs, kwargs...), levels_elicitor)

    # Create SLF
    SLF(model, probes, mean)
end

# Define `smooth` function that a takes SLF as argument with default params
# function Utilities.smooth(::SLF, x, y)
#     # Guestimate reasonable starting places for parameter values

#     # Smooth result and return x̂ and ŷ
#     smooth(x, y, logistic; x₀=[-7.5, 0.15, -70.0, 200.0])
# end

function Utilities.interpolate(::SLF, x, y)
    # Interpolate the response using a linear interpolation
    itp = linear_interpolation(x, y, extrapolation_bc=Line())
    x̂ = LinRange(extrema(x)..., 1000)
    ŷ = itp(x̂)
    return x̂, ŷ
end

# Define `level` to point to elicitor levels
Utilities.level(sim::SLF) = map(x -> level(contra(x)), sim.probes) 
probelevel(sim::SLF) = level(ipsi(sim.probes[1]))  # new function `probe` to point to probe level

# Overload viz
function Utilities.viz(
    sim::SLF; 
    config=Config(), 
    fig=Figure(), 
    ax=Axis(fig[1, 1]), 
    color=:black,
    linestyle=:solid,
    marker=:circle,
    interp=true,
    kwargs...
)
    μ = @memo config simulate(sim)
    scatter!(ax, level(sim), μ; color=color, marker=marker)
    if interp
        x̂, ŷ = Utilities.interpolate(sim, level(sim), μ)
        lines!(ax, x̂, ŷ; color=:red, linestyle=linestyle)
    else
        lines!(ax, level(sim), μ; color=color, linestyle=linestyle)
    end 
    vlines!(ax, threshold(sim; config=config, interp=interp); color=color, linestyle=:dash)
    xlims!(ax, 0.0, 80.0)
    ylims!(ax, 0.0, 280.0)
    display(fig)
end

"""
    threshold(sim::SLF[; config=Config()])

Compute threshold for deterministic SLF as level achieving reduction equal to 10% of maximum rate
"""
function Utilities.threshold(sim::SLF; config=Config(), interp=false)
    μ = @memo config simulate(sim) 
    if interp
        l, μ = Utilities.interpolate(sim, level(sim), μ)
    else
        l = level(sim)
    end
    supp_threshold(l, μ)
end

function supp_threshold(l::AbstractVector, μ::AbstractVector)
    # Identify indices of levels satisfying criterion
    idxs = findall(μ .< 0.95*maximum(μ))

    # If we didn't find any, return NaN, otherwise return level at first idx
    if isempty(idxs)
        return NaN
    else
        return l[first(idxs)]
    end
end

"""
    magnitude(sim::SLF[; config=Config()])

Compute magnitude for SLF as % reduction in rate comparing lowest to highest levels tested
"""
function magnitude(sim::SLF; config=Config())
    # Compute responses
    μ = @memo config simulate(sim)

    # Identify lowest and highest levels
    idx_low = argmin(level(sim))
    idx_high = argmax(level(sim))
  
    # Return change in rate in percentage
    return 100 * (1 - μ[idx_high] / μ[idx_low])
end
