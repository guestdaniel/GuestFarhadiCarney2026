# Export custom CRLB types
export CRLB, CRLBPop, CRLBAI, CRLBAI_Detection, threshold

# Declare abstract type to cover CRLB simulations
abstract type AbstractCRLB <: Simulation end

# Function to convert the output of simulate(crlb) to threshold
function threshold(sim::AbstractCRLB, args...)
    sqrt(simulate(sim, args...))
end

"""
    CRLB <: Simulation

Computes time-averaged Cramer-Rao lower bound with respect to specified parameter.

Note that this CRLB code is designed to work with models that have a single readout channel,
i.e., `ismultichannel(model) == false` because `length(model.coi) == 1`. Note currently
that time is ignored in this calculation.
"""
struct CRLB{M <: Model, S <: AbstractStimulus} <: AbstractCRLB
    model::M
    stim1::S
    stim2::S
    parameterfunc::Function
    summaryfunc::Function
    σ::Float64
    N::Float64
end

function simulate(sim::CRLB)
    # Simulate both responses
    r1 = sim.summaryfunc(compute(sim.model, sim.stim1)[idxswin(sim.stim1)])
    r2 = sim.summaryfunc(compute(sim.model, sim.stim2)[idxswin(sim.stim2)])

    # Determine delta from stimuli objects
    δ = sim.parameterfunc(sim.stim2) - sim.parameterfunc(sim.stim1)

    # Based on the stage of the model, send to either `crlb_poisson` or `crlb_normal`
    # `crlb_normal` is a special case used to analyze BM and IHC responses, where they can
    # be treated as normal distributions after log power transformation. In this case,
    # sim.σ is used as the internal noise parameter that limits the precision of the 
    # estimator.
    if islogout(sim.model)
        crlb_normal(r1, r2, δ, sim.σ, sim.N)
    else
        crlb_poisson(r1, r2, δ, sim.N)
    end
end

"""
    CRLBAI <: Simulation

Computes all-information Cramer-Rao lower bound with respect to specified parameter.

Note that this CRLB code is designed to work with models that have a single readout channel,
i.e., `ismultichannel(model) == false` because `length(model.coi) == 1`. 
"""
struct CRLBAI{M <: Model, S <: AbstractStimulus} <: AbstractCRLB
    model::M
    stim1::S
    stim2::S
    parameterfunc::Function
    σ::Float64
    N::Float64
end

function simulate(sim::CRLBAI)
    # Simulate both responses
    r1 = compute(sim.model, sim.stim1)[idxswin(sim.stim1)]
    r2 = compute(sim.model, sim.stim2)[idxswin(sim.stim2)]

    # Determine delta from stimuli objects
    δ = sim.parameterfunc(sim.stim2) - sim.parameterfunc(sim.stim1)

    # Based on the stage of the model, send to either `crlb_poisson` or `crlb_normal`
    # `crlb_normal` is a special case used to analyze BM and IHC responses, where they can
    # be treated as normal distributions after log power transformation. In this case,
    # sim.σ is used as the internal noise parameter that limits the precision of the 
    # estimator.
    if islogout(sim.model)
        error("All-information CRLB not implemented for normal distributions")
    else
        crlb_ai_poisson(r1, r2, δ, sim.N, samprate(sim.model))
    end
end

"""
    CRLBAI_Detection <: Simulation

Computes all-information optimal detection threshold with respect to specified parameter.

Note that this CRLB code is designed to work with models that have a single readout channel,
i.e., `ismultichannel(model) == false` because `length(model.coi) == 1`. 
"""
struct CRLBAI_Detection{M <: Model, S <: AbstractStimulus} <: AbstractCRLB
    model::M
    stim_standard::S
    stims_target::Vector{S}
    parameterfunc::Function
    σ::Float64
    N::Float64
end

function _simulate(sim::CRLBAI_Detection)
    # Loop over pairs of standard stimulus and target stimuli and get sensitivity for each
    Q = map(sim.stims_target) do stim_target
        # Simulate both responses
        r1 = compute(sim.model, sim.stim_standard)[idxswin(sim.stim_standard)]
        r2 = compute(sim.model, stim_target)[idxswin(stim_target)]

        # Based on the stage of the model, send to either `crlb_poisson` or `crlb_normal`
        # `crlb_normal` is a special case used to analyze BM and IHC responses, where they can
        # be treated as normal distributions after log power transformation. In this case,
        # sim.σ is used as the internal noise parameter that limits the precision of the 
        # estimator.
        if islogout(sim.model)
            error("All-information CRLB not implemented for normal distributions")
        else
            optimal_sensitivity_poisson(r1, r2, sim.N, samprate(sim.model))
        end
    end
end

function simulate(sim::CRLBAI_Detection)
    Q = _simulate(sim)
    threshold(sim, axis(sim), Q)
end

axis(sim::CRLBAI_Detection) = map(x -> sim.parameterfunc(x), sim.stims_target)

function threshold(sim::CRLBAI_Detection, x, y)
    idx = findfirst(>(1), y)
    if isnothing(idx)
        return NaN  # No threshold found
    else
        return x[idx]  # Return the x value at the first index where y exceeds 1
    end
end

"""
    CRLBPop <: Simulation

Computes time-averaged Cramer-Rao lower bound with respect to specified parameter in population response.

Note that this CRLB code is designed to work with models that have multiple readout channels,
i.e., `ismultichannel(model) == true` because `length(model.coi) > 1`. Note currently
that time is ignored in this calculation.
"""
struct CRLBPop{M <: Model, S <: AbstractStimulus} <: AbstractCRLB
    model::M
    stim1::S
    stim2::S
    parameterfunc::Function
    summaryfunc::Function
    σ::Float64
    N::Float64
end

function simulate(sim::CRLBPop)
    # Simulate both responses (will be n_chan-length vectors)
    r1 = compute(sim.model, sim.stim1)
    r2 = compute(sim.model, sim.stim2)
    r1 = map(r1) do chan 
        sim.summaryfunc(chan[idxswin(sim.stim1)])
    end
    r2 = map(r2) do chan 
        sim.summaryfunc(chan[idxswin(sim.stim2)])
    end

    # Determine delta from stimuli objects
    δ = sim.parameterfunc(sim.stim2) - sim.parameterfunc(sim.stim1)

    # Based on the stage of the model, send to either `crlb_poisson` or `crlb_normal`
    if islogout(sim.model)
        crlb_normal(r1, r2, δ, sim.σ, sim.N)
    else
        crlb_poisson(r1, r2, δ, sim.N)
    end
end

"""
    crlb_poisson(λ1::Float, λ2::Float, Δθ::Float, N::Float, T::Float)
    crlb_poisson(λ1::Vector, λ2::Vector, Δθ::Float, N::Float, T::Float)

Compute Cramer-Rao lower bound w.r.t. θ using finite-difference method for Poisson dist.

Assume that λ is the rate-parameter of a Poisson distribution, and that λ is the function of
an unknown parameter θ s.t. λ = f(θ). Supoose further that we know two values of λ, λ(θ) and
λ(θ + ΔΘ) (denoted in the function signature with `λ1` and `λ2` respectively) where Δθ is
small. Assume that T is the interval of time over which we observed this process, and N is the
number of IID processes we observe. We wish to approximate the Cramer-Rao lower bound for
estimating λ, which is:

CRLB(θ) = 1/FI(θ) where FI(θ) is the Fisher information for θ

FI(θ) = N * T * 1/λ(θ) * (dλ(θ)/dθ)^2

This code implements said approximation. If the λ are scalar-valued, then we apply the
equation above. If the λ are vector-valued, we calculate the Fisher information for each
element and then combine via summation before computing the reciprocal. We can do this
because we assume that responses from different processes are independent.
"""
function crlb_poisson(λ1::AbstractFloat, λ2::AbstractFloat, Δθ::AbstractFloat, N::AbstractFloat=1.0, T::AbstractFloat=1.0)
    # Compute partial derivative of λ w.r.t. θ
    ∂θ = (λ2 - λ1) / Δθ

    # Compute Fisher information
    fi = N * T * (1/(λ1 + eps())) * ∂θ^2

    # Return CRLB as reciprocal of Fisher information
    1/fi
end

function crlb_poisson(Λ1::AbstractVector, Λ2::AbstractVector, Δθ::AbstractFloat, N::AbstractFloat=1.0, T::AbstractFloat=1.0)
    # Map through elements of Λ1/Λ2 (implicitly, responses from different channels)
    fi = map(zip(Λ1, Λ2)) do (λ1, λ2)
        # Compute partial derivative of λ w.r.t. θ
        ∂θ = (λ2 - λ1) / Δθ

        # Compute Fisher information 
        fi = N * T * (1/(λ1 + eps())) * ∂θ^2
        return fi
    end
    1/sum(fi)  # combine Fisher information via summation and take reciprocal
end

"""
    crlb_ai_poisson(λ1::Vector, λ2::Vector, Δθ::Float, N::Float)
    crlb_ai_poisson(λ1::Vector{Vector}, λ2::Vector{Vector}, Δθ::Float, N::Float)

Compute all-information Cramer-Rao lower bound w.r.t. θ using finite-difference method for
Poisson dist.

Assume that λ(t) is the rate-parameter of a nonhomogeneous Poisson distribution, and that λ
is the function of an unknown parameter θ s.t. λ(t) = f(θ, t). Supoose further that we
know two values of λ(t), λ(θ, t) and λ(θ + ΔΘ, t) (denoted in the function signature with
`λ1` and `λ2` respectively). We wish to approximate the Cramer-Rao lower bound for
estimating λ, which is the reciprocal of the Fisher information, which is:

FI(θ) = N * ∫_0^T 1/λ(θ, t) [dλ(θ, t)/dθ]^2 dt

This code implements said approximation with finite differencing at each time step. The
value of eps() is added to avoid division by zero. If the inputs are vector-valued, we apply
the above equation, returning the CRLB. If the inputs are vector-of-vector-valued, we repeat
the process for each element and then combine.
"""
function crlb_ai_poisson(λ1::Vector{Float64}, λ2::Vector{Float64}, Δθ::Float64, N::Float64=1.0, fs::Float64=100e3)
    # Compute partial derivative at each time step w.r.t. θ
    ∂θ = @. (λ2 - λ1) / Δθ

    # Compute Fisher information via expression above
    fi_inside = 1.0 ./ (λ1 .+ eps()) .* ∂θ .^ 2
    fi = N * 1/fs * sum(fi_inside)

    # Return CRLB as reciprocal of fisher information
    1/fi
end

function crlb_ai_poisson(Λ1::Vector{Vector{Float64}}, Λ2::Vector{Vector{Float64}}, Δθ::Float64, N::Float64=1.0, fs::Float64=100e3)
    # Map through elements of Λ1/Λ2 (implicitly, responses from different channels)
    fi = map(zip(Λ1, Λ2)) do (λ1, λ2)
        # Compute partial derivative at each time step w.r.t. θ
        ∂θ = @. (λ2 - λ1) / Δθ

        # Compute Fisher information via expression above
        fi_inside = 1.0 ./ (λ1 .+ eps()) .* ∂θ .^ 2
        fi = N * 1/fs * sum(fi_inside)
        return fi
    end

    # Return CRLB as reciprocal of fisher information
    1/sum(fi) 
end

"""
    optimal_sensitivity_poisson_viaclrb(λ1, λ2, Δθ, fs)

Computes the sensitivity for optimal decoding of nonhomogeneous Poisson.

Assume that λ(t, θ) is the time-varying rate parameter of a nonhomogeneous Poisson
distribution. Heinz et al. (2001) show that the squared normalized sensitivity of an optimal
estimator of θ performing discrimination between θ and θ + Δθ is given by:

(δ'_θ)^2 = ∫_0^T 1/λ(t, θ) [∂λ(t, θ)/∂θ]^2 dt 

Note that the inverse square root of this quantity is the threshold of the ideal observer
for discrimination between θ and θ + Δθ, because this is the size of Δθ where sensitivity
is equal to 1. We can transform from squared normalized sensitivity by observing that:

d'^2 = (Δθ)^2 δ'_θ^2  (2.19 of Heinz 2001)

And then taking the square root of both sides gives us:

d' = Δθ δ'_θ

Note that this approximation, depending on use case, may only valid for small Δθ! (However,
it seems to agree very nicely with `optimal_sensitivity_poisson` below, when using a
very low placeholder for -Inf in detection contexts)
"""
function optimal_sensitivity_poisson_viacrlb_ai(λ1::AbstractVector, λ2::AbstractVector, Δθ::AbstractFloat, N::AbstractFloat=1.0, fs::AbstractFloat=100e3)
    ∂θ = @. (λ2 - λ1 + eps()) / Δθ       # derivative with repsect to θ
    fi = @. (N * ∂θ ^ 2) / (λ1 + eps())  # Fisher information for Poisson
    Δθ * sqrt(sum(fi .* 1/fs))           # integrate using rectangular approximation to integral and take square root
end

function optimal_sensitivity_poisson_viacrlb_rp(λ1::AbstractFloat, λ2::AbstractFloat, Δθ::AbstractFloat, N::AbstractFloat=1.0, dur::AbstractFloat=1.0)
    Δθ * sqrt(1/crlb_poisson(λ1, λ2, Δθ, N, dur))
end

"""
    optimal_sensitivity_poisson_rp(λ1::Float, λ2::Float, Δθ, fs)

Computes the sensitivity for optimal decoding of Poisson interval count.

Assume that λ(θ) is the rate parameter of a homogeneous Poisson
distribution. Heinz (2001) shows that the the squared sensitivity of an optimal
estimator of θ, Y, performing discrimination between θ1 and θ2 is given by:
    Q(Y) = [E(Y|θ2) - E(Y|θ1)]^2/var(Y)

It is straightforward to derive from first principles that for a Poisson distribution where
Y is the maximum likelihood estimator, that the sensitivity Q is equal to
    Q(Y) = (λ2-λ1)^2/(λ1/n)

For consistency with typical d' estimates, we return the square root of this quantity.

TODO: figure out how time factors into the equation here!
"""
function optimal_sensitivity_poisson_rp(λ1::AbstractFloat, λ2::AbstractFloat, N::AbstractFloat=1.0, T::AbstractFloat=1.0)
    # Calculate numerator
    num = (λ2 - λ1)^2

    # Now calculate denominator
    den = λ1/N

    return sqrt(num / den)
end

function optimal_sensitivity_poisson_rp(Λ1::Vector{Float64}, Λ2::Vector{Float64}, N::AbstractFloat=1.0, T::AbstractFloat=1.0)
    # Map over each channel of Λ1/Λ2 and compute Q
    Q = map(zip(Λ1, Λ2)) do (λ1, λ2)
        # Calculate numerator
        num = (λ2 - λ1)^2

        # Now calculate denominator
        den = λ1/N

        return num/den
    end
    return sqrt(sum(Q))  # combine sensitivities via summation and take square root to turn into d'
end


"""
    optimal_sensitivity_poisson_ai(λ1::Vector, λ2::Vector, Δθ, fs)

Computes the sensitivity for optimal decoding of nonhomogeneous Poisson data.

Assume that λ(t, θ) is the time-varying rate parameter of a nonhomogeneous Poisson
distribution. Heinz (2001) shows that the the squared sensitivity of an optimal
estimator of θ, Y, performing discrimination between θ1 and θ2 is given by:
    Q(Y) = [E(Y|θ2) - E(Y|θ1)]^2/var(Y)

As explained in Eq. 5.22 of Heinz (2001), it turns out that:
    E(Y) = ∫_0^T ln(λ(t|θ2)/λ(t|θ1)) λ(t) dt
    Var(Y) = ∫_0^T [ln(λ(t|θ2)/λ(t|θ1))]^2 λ(t) dt

Q then simplifies to:
    (∫_0^T ln(λ(t|θ2)/λ(t|θ1)) [λ(t|θ2) - λ(t|θ1)] dt)^2 / (∫_0^T [ln(λ(t|θ2)/λ(t|θ1))]^2 λ(t) dt)

We return, for consistency, the square root of this quantity, which is the sensitivity of
the ideal observer Y.
"""
function optimal_sensitivity_poisson_ai(λ1::Vector{Float64}, λ2::Vector{Float64}, N::AbstractFloat=1.0, fs::AbstractFloat=100e3, o::AbstractFloat=1.0)
    # First, calculate the quantity ln(λ(t|θ2)/λ(t|θ1))
    g = @. log((λ2 + o)/(λ1 + o))  # vector of point-by-point calculation for log ratio of rates

    # Now calculate numerator
    num = (N * sum(@. 1/fs * g * (λ2 - λ1)))^2

    # Now calculate denominator
    den = N * sum(@. 1/fs * g^2 * λ1)

    return sqrt(num / den)
end

function optimal_sensitivity_poisson_ai(Λ1::Vector{Vector{Float64}}, Λ2::Vector{Vector{Float64}}, N::AbstractFloat=1.0, fs::AbstractFloat=100e3, o::AbstractFloat=1.0)
    # Map over each channel of Λ1/Λ2 and compute Q
    Q = map(zip(Λ1, Λ2)) do (λ1, λ2)
        # First, calculate the quantity ln(λ(t|θ2)/λ(t|θ1))
        g = @. log((λ2 + o)/(λ1 + o))  # vector of point-by-point calculation for log ratio of rates

        # Now calculate numerator
        num = (N * sum(@. 1/fs * g * (λ2 - λ1)))^2

        # Now calculate denominator
        den = N * sum(@. 1/fs * g^2 * λ1)
        return num/den
    end

    return sqrt(sum(Q))
end

"""
    crlb_normal(μ1::AbstractFloat, μ2::AbstractFloat, Δθ::AbstractFloat, σ::AbstractFloat)
    crlb_normal(M1::AbstractVector, M2::AbstractVector, Δθ::AbstractFloat, σ::AbstractFloat)

Compute Cramer-Rao lower bound w.r.t. θ using finite-difference method for normal dist.

Assume that μ is the rate-parameter of a normal distribution, and that μ is the function
of an unknown parameter θ s.t. μ = f(θ). Supoose further that we observe two values of 
μ, μ(θ) and μ(θ + ΔΘ) (denoted in the function signature with `μ1` and `μ2` respectively). 
We wish to approximate the Cramer-Rao lower bound for estimating μ, which is:

CRLB(θ) = σ^2/(dμ(θ)/dθ)^2

This code implements said approximation. Note that observation time is currently ignored!
σ is treated as a known parameter in this function and is equivalent to an internal noise
parameter for most applications. If the μ are scalar-valued, then we apply the equation
above. If the μ are vector-valued, we repeat the process for each element and then combine
results via summation before returning.
"""
function crlb_normal(μ1::AbstractFloat, μ2::AbstractFloat, Δθ::AbstractFloat, σ::AbstractFloat, N::Float64=1.0)
    ∂θ = (μ2 - μ1) / Δθ  # derivative with respect to θ
    σ^2 / (N * ∂θ^2)     # reciprocal of Fisher information for normal
end

function crlb_normal(M1::AbstractVector, M2::AbstractVector, Δθ::AbstractFloat, σ::AbstractFloat, N::Float64=1.0)
    fi = map(zip(M1, M2)) do (μ1, μ2)
        ∂θ = (μ2 - μ1) / Δθ  # derivative with respect to θ
        (N * ∂θ^2) / σ^2     # Fisher information for Poisson
    end
    1/sum(fi)  # combine Fisher information via summation and take reciprocal
end

"""
    viz(sims::Vector{::Model, ::AbstractStimulus}; kwargs...)
"""
function viz(
    sims::Vector{<:AbstractCRLB};
    config=Config(),
    type="ldl",
    kwargs...
)
    # First, simulate the result and extract a suitable x-axis
    θ = @memo config threshold(sims)
    x = map(x -> x.parameterfunc(x.stim1), sims)

    # Sent to correct type of plot
    if type == "ldl"
        fig, ax = plot_ldl(x, θ; kwargs...)
    elseif type == "amdl"
        fig, ax = plot_ldl(map(x -> level(x.stim1), sims), θ; kwargs...)
    else
        error("Unknown type: $type")
    end

    display(fig)
    return fig, ax
end

"""
    viz(sim::CRLBAI_Detection; kwargs...)
"""
function viz(
    sim::CRLBAI_Detection;
    config=Config(),
    fig=Figure(),
    ax=Axis(fig[1, 1]; yscale=log10),
    kwargs...
)
    # Simulate d' at every point
    d_prime = @memo config Utilities._simulate(sim)
    scatter!(ax, axis(sim), d_prime; color=:black, kwargs...)
    lines!(ax, axis(sim), d_prime; color=:black, linewidth=2.0, kwargs...)

    # Extract threshold
    θ = threshold(sim, axis(sim), d_prime)
    scatter!(ax, [θ], [1.0]; color=:red)

    display(fig)
    return fig, ax
end

# `viz` method for vectors of CRLB simulations, presumed to each map to a vector of thresholds
# and be plotted against axis
function viz(
    simsvec::Vector{Vector{T}};
    config=Config(),
    type="ldl",
    colors=colorschemes[:Dark2_8],
    fig=Figure(; size=(450, 400)),
    ax=Axis(fig[1, 1]; yscale=type == "ldl" ? log10 : identity),
    kwargs...
) where {T <: AbstractCRLB}
    # First, simulate the result and extract a suitable x-axis
    results = pmap(enumerate(simsvec)) do (idx, sims)
        θ = @memo config threshold(sims)
        x = map(x -> x.parameterfunc(x.stim1), sims)
        x, θ
    end

    # Sent to correct type of plot
    map(enumerate(results)) do (idx, (x, θ))
        if type == "ldl"
            plot_ldl(x, θ; fig=fig, ax=ax, setup_axis=idx==1, color=colors[idx], kwargs...)
        end
    end

    display(fig)
    return fig, ax
end

# `viz` method for vectors of CRLB*_Detection simulations, presumed to each map to a single
# threshold and be plotted against sound level
function viz(
    sims::Vector{T};
    config=Config(),
    fig=Figure(; size=(450, 400)),
    ax=Axis(fig[1, 1]),
    color=:black,
    setup_axis=true,
    linewidth=1.0,
    kwargs...
) where {T <: CRLBAI_Detection}
    # First, simulate the result and extract a suitable x-axis
    results = pmap(sims) do sim
        @memo config simulate(sim)
    end
    x = map(x -> level(x.stim_standard), sims)

    # If setup_axis, add ticks, labels, limits, etc.
    if setup_axis
        # Set ticks
        lmin = floor(minimum(x); digits=-1)
        lmax = ceil(maximum(x); digits=-1)
        ax.xticks = lmin:20:lmax
        ax.xminorticks = lmin:5:lmax
        ax.xminorticksvisible = true

        # Set lims
        xlims!(ax, lmin - 5.0, lmax + 5.0)

        # Set labels
        ax.xlabel = "Level (dB SPL)"
        ax.ylabel = "AMDL (dB)"
    end

    # Plot against each other
    lines!(ax, x, results; color=color, linewidth=linewidth)

    display(fig)
    return fig, ax
end

function viz(
    sims::Vector{Vector{T}};
    config=Config(),
    fig=Figure(; size=(450, 400)),
    colors=colorschemes[:Dark2_8],
    ax=Axis(fig[1, 1]),
    kwargs...
) where {T <: CRLBAI_Detection}
    map(zip(sims, colors)) do (sim, c)
        viz(sim; config=config, fig=fig, ax=ax, color=c)
    end
    display(fig)
    return fig, ax
end

function plot_ldl(
    l, 
    ldl;
    fig=Figure(; size=(450, 350)),
    yscale=log10,
    ax=Axis(fig[1, 1]; yscale=yscale),
    setup_axis=true,
    linestyle=:solid,
    color=:black,
    linewidth=1.0,
    ylims=[],
    pad=1,
)
    # If setup_axis, add ticks, labels, limits, etc.
    if setup_axis
        # Set ticks
        lmin = floor(minimum(l); digits=-1)
        lmax = ceil(maximum(l); digits=-1)
        ax.xticks = lmin:20:lmax
        ax.xminorticks = lmin:5:lmax
        ax.xminorticksvisible = true
        ax.yticks = autolog10ticks(ldl)
        ax.yminorticksvisible = true
        ax.yminorticks = IntervalsBetween(9)

        # Set lims
        xlims!(ax, lmin - 5.0, lmax + 5.0)
        if isempty(ylims)
            ylims!(ax, autolog10lims(ldl; pad=pad)...)
        else
            ylims!(ax, ylims)
        end

        # Set labels
        ax.xlabel = "Level (dB SPL)"
        ax.ylabel = "LDL (dB)"
    end

    # Plot data with color
    lines!(ax, l, ldl; color=color, linestyle=linestyle, linewidth=linewidth)

    display(fig)
    return fig, ax
end

function plot_amdl(
    l,
    amdl;
    fig=Figure(; size=(450, 350)),
    ax=Axis(fig[1, 1]),
    setup_axis=true,
    linestyle=:solid,
    color=:black,
    linewidth=1.0,
    ylims=[],
    pad=NaN,
)
    # If setup_axis, add ticks, labels, limits, etc.
    if setup_axis
        # Set ticks
        lmin = floor(minimum(l); digits=-1)
        lmax = ceil(maximum(l); digits=-1)
        ax.xticks = lmin:20:lmax
        ax.xminorticks = lmin:5:lmax
        ax.xminorticksvisible = true

        # Set lims
        xlims!(ax, lmin - 5.0, lmax + 5.0)

        # Set labels
        ax.xlabel = "Level (dB SPL)"
        ax.ylabel = "AMDL (dB)"
    end

    # Plot data
    lines!(ax, l, amdl; color=color, linestyle=linestyle, linewidth=linewidth)
    if !isempty(ylims)
        ylims!(ax, ylims...)
    end

    display(fig)
    return fig, ax
end