export PFObserver, PFTemplateObserver, summarize, fit, mahalanobis, logistic, invlogistic, logistic_fit, fit

abstract type PF <: Simulation end

"""
    PFObserver

Type for estimating psychometric function from sequence of excitation patterns
"""
@with_kw struct PFObserver{S, M, N} <: PF where {S <: AbstractStimulus, M <: Model, N <: Int}
    patterns::Vector{NTuple{N, AvgPattern{S, M}}}
    model::M=patterns[1][1].model
    preprocessor::Function=pre_nothing
    observer::Function=obs_maxrate
    n_step::Int64=length(patterns)
    n_trial::Int64=patterns[1][1].n_rep
    n_interval::Int64=2  # todo fix this
    θ::Symbol
    θ_low=getfield(patterns[1][end].stimuli[1], θ)
    θ_high=getfield(patterns[end][end].stimuli[1], θ)
end

function simulate(s::PFObserver, config::Config=Default())
    # First, we go through each step of the psychometric function and simulate responses and
    # apply the observer to make a decision about each response pair. Recall that each 
    # "pattern" is one step along the psychometric function abscissa encoding information
    # about how to simulate each trial
    map(enumerate(s.patterns)) do (idx, pattern_tuple)
        # Print info
        @info "Evaluating PF, step $idx of $(s.n_step)"

        # Each element of pattern is itself an AvgPattern that needs to be evaluated or 
        # loaded from disk. There should be N AvgPatterns, where N is the number of 
        # intervals, and each AvgPattern should have n_trial elements
        outs = map(pattern_tuple) do p
            @memo config simulate(p)
        end

        # We can apply a preproessing function, stored as field s.preprocessor, if we 
        # want to do any normalization, artifact rejection, etc.
        outs = s.preprocessor(outs)

        # Loop through each trial, permute order, and make decision based on observer...
        # Then, return correct/incorrect for each trial. μs is a vector of length N 
        # containing the response to each interval on the current iterated trial
        map(zip(outs...)) do μs
            order = shuffle(1:s.n_interval)
            decision = s.observer(μs[order]...)
            order[decision] == s.n_interval ? 1 : 0  # convert to correct/incorrect
        end
    end
end

"""
    PFTemplateObserver

Type for estimating psychometric function from template and sequence of excitation patterns
"""
@with_kw struct PFTemplateObserver{S, M, N} <: PF where {S <: AbstractStimulus, M <: Model, N <: Int}
    template::AvgPattern{S, M}
    patterns::Vector{NTuple{N, AvgPattern{S, M}}}
    model::M=patterns[1][1].model
    observer::Function=obs_mahalanobis
    n_step::Int64=length(patterns)
    n_rep_template::Int64=template.n_rep
    n_trial::Int64=patterns[1][1].n_rep
    n_interval::Int64=2  # todo fix this
    θ::Symbol
    θ_low=getfield(patterns[1][end].stimuli[1], θ)
    θ_high=getfield(patterns[end][end].stimuli[1], θ)
end

function simulate(s::PFTemplateObserver, config::Config=Default())
    # First, simulate the template responses
    @info "Loading or estimating template for PF"
    template = @memo config simulate(s.template)
    μ = mean(template)
    Σ = cov(template)

    # Next, we go through each step of the psychometric function and simulate responses and
    # apply the observer to make a decision about each response pair. Recall that each 
    # "pattern" is one step along the psychometric function abscissa encoding information
    # about how to simulate each trial
    map(enumerate(s.patterns)) do (idx, pattern_tuple)
        # Print info
        @info "Evaluating PF, step $idx of $(s.n_step)"

        # Each element of pattern is itself an AvgPattern that needs to be evaluated or 
        # loaded from disk. There should be N AvgPatterns, where N is the number of 
        # intervals, and each AvgPattern should have n_trial elements
        outs = map(pattern_tuple) do p
            @memo config simulate(p)
        end

        # Loop through each trial, permute order, and make decision based on observer...
        # Then, return correct/incorrect for each trial. μs is a vector of length N 
        # containing the response to each interval on the current iterated trial
        map(zip(outs...)) do μs
            order = shuffle(1:s.n_interval)
            decision = s.observer(μ, Σ, μs[order]...)
            order[decision] == s.n_interval ? 1 : 0  # convert to correct/incorrect
        end
    end
end

function mahalanobis(μ::Vector{Float64}, Σ::Matrix{Float64}, x::Vector{Float64})
    sqrt((μ .- x)' * inv(Σ) * (μ .- x))
end

function mahalanobis(template::CovMatrix, x::Vector{Float64})
    μ = mean(template)
    Σ = cov(template)
    sqrt((μ .- x)' * inv(Σ) * (μ .- x))
end

function obs_mahalanobis(μ::Vector{Float64}, Σ::Matrix{Float64}, x::Vararg{Vector{Float64}})
    argmax(map(_x -> mahalanobis(μ, Σ, _x), x))
end

function obs_maxrate(x::Vararg{Vector{Float64}})
    argmax(map(mean, x))
end

pre_nothing(args) = args

# Generic function for implementing ID for a psychoemtric function
function id(p::PF; accesses=nothing, connector="_", kwargs...)
    # Get id components corresponding to interpretable parts (model, n_rep, summarizer)
    main_part = savename(
        string(typeof(p)),
        p; 
        accesses=accesses === nothing ? fieldnames(typeof(p)) : accesses,
        allowedtypes=(
            Real, 
            String, 
            Symbol, 
            Function,
            Component,
            Audiogram,
        ), 
        connector=connector,
        kwargs...
    )

    # Get id component for first stimulus in stimuli
    stim_1 = id(p.patterns[1][end].stimuli[1])
    stim_2 = id(p.patterns[end][end].stimuli[1])

    # Return combination
    return main_part * connector * stim_1 * connector * stim_2
end

# Functions to handle fitting logistic functions to psychometric function data
"""
    logistic(x, λ; L, offset)

Returns the logistic function value evaluated at X

# Arguments
- `x`: Argument
- `λ`: Midpoint and slope of the logistic function (in a tuple/array)
- `L=1.0`: Maximum value of the logistic function (before offset)
- `offset=0.0`: Offset of the logistic function
"""
function logistic(x, λ; L=1.0, offset=0.0)
    L/(1 + exp(-λ[2]*(x-λ[1]))) + offset
end

function invlogistic(y, λ; L=1.0, offset=0.0)
    log(L/(y-offset) - 1)/(-λ[2]) + λ[1]
end

"""
    logistic(x, t, s; L, offset)

Returns the logistic function value evaluated at X

# Arguments
- `x`: Argument
- `t`: Midpoint of the logistic function
- `s`: Slope of the logistic function
- `L=1.0`: Maximum value of the logistic function (before offset)
- `offset=0.0`: Offset of the logistic function
"""
function logistic(x, t, s; L=1.0, offset=0.0)
    L/(1 + exp(-s*(x-t))) + offset
end

# Vector-valued function that maps from vector x and pair p (threshold and
# slope) to logistic output
@. logistic_fit(x, p) = logistic(x, p[1], p[2]; L=0.5, offset=0.5)

# Vector-valued function that maps from vector x, vector p1 (threshold), and
# vector p2 (slope) to logistic output
logistic_predict(x, p1, p2) = logistic.(x, p1, p2; L=0.5, offset=0.5)

# Function that handles fitting using above functions automatically
function fit_psychometric_function(x, y)
    curve_fit(
        logistic_fit,
        x,
        y,
        [0.0, 1.0];
        lower=[-40.0, 0.01],
        upper=[20.0, 1.0]
    )
end

function fit(::PF, x, y)
    fit_psychometric_function(x, map(mean, y))
end