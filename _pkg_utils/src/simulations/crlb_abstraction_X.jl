# This source file implements the abstraction layer for Cramer-Rao lower bound simulations
# based on detection and/or discrimination tasks. The basic philosophy of the abstraction 
# layer is as follows:
# - A Task is type that encodes information about...
# - A Trial is a type that can wrap around a Task to produce...
# - An Observer is a type that defines how a Task is to be performed and can be interrogated
#   with simulation functions like `threshold`

# Export types
export TaskObserverFeature

export Discrimination, Detection, DetectionAdaptive

export ScalarFeature, VectorFeature

export PoissonIdealObserver, NormalIdealObserver, ConstantStimulusObserver

# Define AbstractTask type
abstract type AbstractTask <: Component end
abstract type AbstractDiscriminationTask <: AbstractTask end
abstract type AbstractDetectionTask <: AbstractTask end

# `trial` returns a standard Stimulus and a target Stimulus for a given parameter Δ
function trial(task::AbstractDiscriminationTask, Δ::Float64)
    task.standard, task.targetfunc(Δ)
end

struct Discrimination{S} <: AbstractDiscriminationTask where {S<:AbstractStimulus}
    standard::S
    targetfunc::Function
end

struct Detection{S} <: AbstractDetectionTask where {S<:AbstractStimulus}
    standard::S
    targets::Vector{S}
    abscissa::Vector{Float64}
end

struct DetectionAdaptive{S} <: AbstractDetectionTask where {S<:AbstractStimulus}
    standard::S
    targetfunc::Function
    pinit::Float64
    pmax::Float64
    pmin::Float64
    step_size_init::Float64
    step_size_factor::Float64
    n_step::Int64
    sign::Float64                  # if sign == -1, then a downward step corresponds to a decrease in the parameter value
end

# Define AbstractFeature type
abstract type AbstractFeature <: Component end
abstract type AbstractVectorFeature <: AbstractFeature end
abstract type AbstractScalarFeature <: AbstractFeature end

# ScalarFeature applies a `featurefunc(r)` that returns a scalar value after subsetting 
# the vector-valued response waveform to the indices specified by `idxsfunc(stim)` and
# windowing it with `windowfunc(length(stm))`.
@with_kw struct ScalarFeature <: AbstractScalarFeature
    idxsfunc::Function = idxswin
    windowfunc::Function = rect
    featurefunc::Function = mean
end

# VectorFeature applies a `transformfunc.(r)` that returns a vector after subsetting 
# the vector-valued response waveform to the indices specified by `idxsfunc(stim)` 
@with_kw struct VectorFeature <: AbstractVectorFeature
    idxsfunc::Function = idxswin
    transformfunc::Function = identity
end

# Define a new method `compute` that takes a model, a stimulus, and a `feature`
function compute(model::Model, stim::AbstractStimulus, feature::AbstractFeature)
    # Simulate response (response could be single vector or vector of vectors in most cases)
    r = compute(model, stim)

    # Branch based on the type of r
    _extract_feature(r, stim, feature)
end

# Method for `_extract_feature` that covers the use case of a ScalarFeature
function _extract_feature(r::Vector{Float64}, stim::AbstractStimulus, feature::AbstractScalarFeature)
    # Subset r to only indices we care about
    r_subset = r[feature.idxsfunc(stim)]

    # Apply window to r_subset
    r_subset = r_subset .* feature.windowfunc(length(r_subset))

    # Extract feature from r_susbet
    feature.featurefunc(r_subset)
end

# Method for `_extract_feature` that covers the use case of a VectorFeature
function _extract_feature(r::Vector{Float64}, stim::AbstractStimulus, feature::AbstractVectorFeature)
    # Subset r to only indices we care about
    r_subset = r[feature.idxsfunc(stim)]

    # Extract feature from r_susbet
    feature.transformfunc.(r_subset)
end

function _extract_feature(r::Vector{Vector{Float64}}, stim::AbstractStimulus, feature::AbstractFeature)
    # For a vector of vector input, map over each element and apply _extract_feature
    map(r) do _r
        _extract_feature(_r, stim, feature)
    end
end

# Define AbstractObserver type
abstract type AbstractObserver <: Component end
abstract type AbstractIdealObserver <: AbstractObserver end

# Define type that combines together a Task + Observer + Feature for convenience
struct TaskObserverFeature{T,O,F} <: Component where {T<:AbstractTask,O<:AbstractObserver,F<:AbstractFeature}
    task::T
    observer::O
    feature::F
end

# Implement flavors of CRLB observers
struct PoissonIdealObserver{M} <: AbstractIdealObserver where {M<:Model}
    model::M
    Δ::Float64
    N::Float64
end

# Calculate threshold for ideal observer with scalar-valued output feature for each channel
function threshold(obs::PoissonIdealObserver, task::AbstractDiscriminationTask, feature::ScalarFeature)
    # Fetch trial from discrimination task object
    standard, target = trial(task, obs.Δ)

    # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
    r1 = compute(obs.model, standard, feature)
    r2 = compute(obs.model, target, feature)

    # Use Cramer-Rao lower bound analysis to determine thresholds (because r1 and r2 are 
    # scalar valued, we employ RP or time-averaged CRLB analysis)
    sqrt(crlb_poisson(r1, r2, obs.Δ, obs.N, dur(standard)))
end

# Calculate threshold for ideal observer with *vector*-valued output feature for each channel
function threshold(obs::PoissonIdealObserver, task::AbstractDiscriminationTask, feature::VectorFeature)
    # Fetch trial from discrimination task object
    standard, target = trial(task, obs.Δ)

    # Simulate both responses (r1 and r2 will be vector-valued, or vector-of-vector-valued)
    r1 = compute(obs.model, standard, feature)
    r2 = compute(obs.model, target, feature)

    # Use Cramer-Rao lower bound analysis to determine thresholds (because r1 and r2 are 
    # vector-valued, we can employ AI CRLB analysis)
    sqrt(crlb_ai_poisson(r1, r2, obs.Δ, obs.N, samprate(obs.model)))
end

# Calculate threshold for ideal-observer with *scalar* valued output feature for each channel
# In this case, we consider the case of a Detection task, versus a Discrimination task, which
# does introduce some complexities. These methods in particular are for the base Detection
# type, which is essentially a constant-stimulus detection task. Sensitivity is estimated
# for each value specified in `task.abscissa` and then threshold is returned as the value
# closest to sensitivity=1.
function _simulate(obs::PoissonIdealObserver, task::Detection, feature::ScalarFeature)
    # Loop over pairs of standard stimulus and target stimuli and get sensitivity for each
    map(task.targets) do target
        # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
        r1 = compute(obs.model, task.standard, feature)
        r2 = compute(obs.model, target, feature)

        optimal_sensitivity_poisson_rp(r1, r2, obs.N, dur(target))
    end
end

# Calculate threshold for ideal-observer with *vector* valued output feature for each channel
# In this case, we consider the case of a Detection task, versus a Discrimination task, which
# does introduce some complexities. These methods in particular are for the base Detection
# type, which is essentially a constant-stimulus detection task. Sensitivity is estimated
# for each value specified in `task.abscissa` and then threshold is returned as the value
# closest to sensitivity=1.
function _simulate(obs::PoissonIdealObserver, task::Detection, feature::VectorFeature)
    # Loop over pairs of standard stimulus and target stimuli and get sensitivity for each
    map(task.targets) do target
        # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
        r1 = compute(obs.model, task.standard, feature)
        r2 = compute(obs.model, target, feature)

        optimal_sensitivity_poisson_ai(r1, r2, obs.N, samprate(obs.model))
    end
end

# Generic method for turning sensitivity Q from `_simulate` into threshold
function threshold(obs::PoissonIdealObserver, task::Detection, feature::AbstractFeature)
    Q = _simulate(obs, task, feature)
    threshold(task, task.abscissa, Q)
end

# Calculate threshold for ideal-observer with *scalar* valued output feature for each channel
# This method extends the logic above for an adaptive sensitivity task, wherein sensitivity
# is estimated between a min and max parameter value using a simple staircase procedure.
# The sensitivity values are then interpolated as before.
function _simulate(obs::PoissonIdealObserver, task::DetectionAdaptive, feature::ScalarFeature)
    # Create variable to store tracked sensitivity and parameter values
    p = Float64[]
    q = Float64[]
    push!(p, task.pinit)

    # Create state variables for loop
    stop = false
    idx = 1
    step = task.step_size_init

    # Loop over pairs of standard stimulus and target stimuli and track sensitivity
    while !stop
        # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
        r1 = compute(obs.model, task.standard, feature)
        r2 = compute(obs.model, task.targetfunc(p[idx]), feature)

        # Compute sensitivity
        push!(q, optimal_sensitivity_poisson_rp(r1, r2, obs.N, samprate(obs.model)))

        # If sensitivity is above 1, we should add sign*step to p[end] to get the new 
        # parameter value; if the sensitivity is below 1, we should add -sign*step to p[end]
        # to get the new parameter value
        if q[end] < 1.0
            push!(p, -task.sign * step + p[end])
        else
            push!(p, task.sign * step + p[end])
        end

        # Now, handle step size — if we have finished all our step sizes, we return!
        step *= task.step_size_factor

        # Increment state variables
        idx += 1

        # Check termination conditions
        if idx > task.n_step
            stop = true
        end

        if (p[end] > task.pmax) || (p[end] < task.pmin)
            stop = true
        end
    end

    # Return p and q
    return p[1:(end-1)], q
end

# Calculate threshold for ideal-observer with *vector* valued output feature for each channel
# This method extends the logic above for an adaptive sensitivity task, wherein sensitivity
# is estimated between a min and max parameter value using a simple staircase procedure.
# The sensitivity values are then interpolated as before.
function _simulate(obs::PoissonIdealObserver, task::DetectionAdaptive, feature::VectorFeature)
    # Create variable to store tracked sensitivity and parameter values
    p = Float64[]
    q = Float64[]
    push!(p, task.pinit)

    # Create state variables for loop
    stop = false
    idx = 1
    step = task.step_size_init

    # Loop over pairs of standard stimulus and target stimuli and track sensitivity
    while !stop
        # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
        r1 = compute(obs.model, task.standard, feature)
        r2 = compute(obs.model, task.targetfunc(p[idx]), feature)

        # Compute sensitivity
        push!(q, optimal_sensitivity_poisson_ai(r1, r2, obs.N, samprate(obs.model)))

        # If sensitivity is above 1, we should add sign*step to p[end] to get the new 
        # parameter value; if the sensitivity is below 1, we should add -sign*step to p[end]
        # to get the new parameter value
        if q[end] < 1.0
            push!(p, -task.sign * step + p[end])
        else
            push!(p, task.sign * step + p[end])
        end

        # Now, handle step size — if we have finished all our step sizes, we return!
        step *= task.step_size_factor

        # Increment state variables
        idx += 1

        # Check termination conditions
        if idx > task.n_step
            stop = true
        end

        if (p[end] > task.pmax) || (p[end] < task.pmin)
            stop = true
        end
    end

    # Return p and q
    return p[1:(end-1)], q
end

function threshold(::AbstractDetectionTask, x, y)
    idx = findfirst(>(1), y)
    if isnothing(idx)
        return NaN  # No threshold found
    else
        return x[idx]  # Return the x value at the first index where y exceeds 1
    end
end

function threshold(obs::PoissonIdealObserver, task::DetectionAdaptive, feature::AbstractFeature)
    # Simulate sensitivity `q` at a range of parameter values `p`
    p, q = _simulate(obs, task, feature)

    # Interpolate using a cubic spline to 1000 values from min to max
    itp = loess(p, q)
    p̂ = LinRange(minimum(p), maximum(p), 1000)
    threshold(task, p̂, predict(itp, p̂))
end

# Extend method `threshold` for TaskObserverFeatures
function threshold(tof::TaskObserverFeature)
    threshold(tof.observer, tof.task, tof.feature)
end

function compute(tof::TaskObserverFeature)
    compute(tof.observer, tof.task, tof.feature)
end

# Implement flavors of CRLB observers
struct NormalIdealObserver{M} <: AbstractIdealObserver where {M<:Model}
    model::M
    Δ::Float64
    σ::Float64
    N::Float64
end

# Calculate threshold for Normal ideal observer with scalar-valued output feature for each channel
function threshold(obs::NormalIdealObserver, task::AbstractDiscriminationTask, feature::ScalarFeature)
    # Fetch trial from discrimination task object
    standard, target = trial(task, obs.Δ)

    # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
    r1 = compute(obs.model, standard, feature)
    r2 = compute(obs.model, target, feature)

    # Use Cramer-Rao lower bound analysis to determine thresholds (because r1 and r2 are 
    # scalar valued, we employ RP or time-averaged CRLB analysis)
    sqrt(crlb_normal(r1, r2, obs.Δ, obs.σ, obs.N))
end

# Convenience plot for showing intermediate results of single TOF simulation (empirical)
function viz(tof::TaskObserverFeature; config=Config(), xscale=log10)
    # Extract abscissa and sensitivity
    x = tof.observer.Δs
    y = @memo config compute(tof)
    th = threshold(tof.observer, tof.task, tof.feature, y)

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=xscale)
    lines!(ax, x, y; color=:black)
    vlines!(ax, [th]; color=:red)

    # Dipslay and return figure
    display(fig)
    return fig, ax
end

# Convenience plots for multiple TOF simulations against an x-axis (default `level`)
function viz(
    sims::Vector{<:TaskObserverFeature},
    parameterfunc=level;
    config=Config(),
    type=nothing,
    kwargs...
)
    # First, simulate the result and extract a suitable x-axis
    θ = pmap(sims) do sim
        @memo config threshold(sim)
    end
    x = map(x -> parameterfunc(x.standard), getfield.(sims, :task))

    # If type is nothing, we can try to guess and auto-select
    if isnothing(type)
        if hasfield(typeof(sims[1].task.standard), :fm)
            type = "amdl"
        else
            type = "ldl"
        end
    end

    # Sent to correct type of plot
    if type == "ldl"
        fig, ax = plot_ldl(x, θ; kwargs...)
    elseif type == "amdl"
        fig, ax = plot_amdl(x, θ; kwargs...)
    else
        error("Unknown type: $type")
    end

    display(fig)
    return fig, ax
end

function viz(
    sims::Vector{Vector{T}},
    parameterfunc=level;
    type="ldl",
    fig=Figure(; size=(450, 400), yscale=type == "ldl" ? log10 : identity),
    ax=Axis(fig[1, 1]),
    colors=colorschemes[:Dark2_8],
    kwargs...
) where {T<:TaskObserverFeature}
    map(enumerate(zip(sims, colors))) do (idx, (_sims, c))
        viz(_sims, parameterfunc; fig=fig, ax=ax, color=c, setup_axis=idx == 1, kwargs...)
    end

    display(fig)
    return fig, ax
end

############################################################################################
# Empirical stuff 

# Define observer type specific to constant stimulus with empirical threshold estimation
struct ConstantStimulusObserver{M} <: AbstractObserver where {M<:Model}
    model::M             # model used here
    Δs::Vector{Float64}  # abscissa values, should always include 0 or -Inf as appropriate
    n_rep::Int64         # how many trials to simulate per value
end

# Calculate sensitivity for empirical observer with scalar-valued output feature
function compute(obs::ConstantStimulusObserver, task::AbstractDiscriminationTask, feature::ScalarFeature)
    # For each level of the abscissa, characterize sensitivity
    map(obs.Δs) do Δ
        # For each rep, compute responses
        # If the model is a single-channel model, then we expect r1 and r2 each to be scalar
        # and therefore targets and standards below will be Vector{<:Float64}
        # If the model is multichannel, however, we need to handle things differently; each
        # response will be a vector of scalars, and therefore targets and standards below will
        # be Vector{Vector{<:Float64}}
        resps = pmap(1:obs.n_rep) do _
            # Fetch trial from discrimination task object
            standard, target = trial(task, Δ)

            # Simulate both responses (r1 and r2 will be scalar-valued, or vector-of-scalar valued)
            r1 = compute(obs.model, standard, feature)
            r2 = compute(obs.model, target, feature)
            return r1, r2
        end

        # Separate resps into targets (second index) and standards (first index)
        targets = getindex.(resps, 2)
        standards = getindex.(resps, 1)

        # Branch based on model type
        if issinglechannel(obs.model)
            return d′(standards, targets)
        else
            targets = [[v[i] for v in targets] for i in 1:length(targets[1])]  # turn into vec{vec{Float64}} where outer vector is CF channels
            standards = [[v[i] for v in standards] for i in 1:length(standards[1])]
            return d′_pop_unbiased(standards, targets)
        end
    end
end

function threshold(obs::ConstantStimulusObserver, task::AbstractDiscriminationTask, feature::ScalarFeature, sens::Vector{Float64})
    # Interpolate sensitivity vs abscissa to find threshold at d' = 1
    a = obs.Δs
    b = sens
    itp = loess(a, b)
    â = LinRange(minimum(a), maximum(a), 1000)
    b̂ = predict(itp, â)
    idx = findfirst(>(1), b̂)
    if idx === nothing
        return NaN
    else
        return â[idx]
    end
end

function threshold(obs::ConstantStimulusObserver, task::AbstractDiscriminationTask, feature::ScalarFeature)
    # Compute sensitivity
    sens = compute(obs, task, feature)
    threshold(obs, task, feature, sens)
end