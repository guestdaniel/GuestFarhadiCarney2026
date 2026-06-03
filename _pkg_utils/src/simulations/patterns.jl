export Pattern, EPattern, DPattern, AvgPattern, FullPattern, DeltaPattern

# Declare abstract types 
abstract type Pattern <: Simulation end
abstract type EPattern <: Pattern end
abstract type DPattern <: Pattern end
#axis(p::Vector{EPattern}, fieldname::Symbol) = map(x -> getfield(x[1], fieldname), p.stimuli)
#axis(p::Vector{DPattern}, fieldname::Symbol) = map(x -> getfield(x[1][1], fieldname), p.stimuli)

# Generic functions
cfaxis(p::Pattern) = p.model.cf

# Define FullPattern (pattern emitting full neurograms for each stimulus)
@with_kw struct FullPattern{S, M} <: EPattern where {S <: AbstractStimulus, M <: Model}
    stimuli::Vector{S}
    model::M
    n_rep::Int64=length(stimuli)
end

function simulate(p::FullPattern)
    @showprogress "Estimating neurogram pattern" pmap(1:p.n_rep) do idx
        compute(p.model, p.stimuli[idx])
    end
end

# Define AvgPattern (pattern emitting summarized neurograms for each stimulus)
@with_kw struct AvgPattern{S, M} <: EPattern where {S <: AbstractStimulus, M <: Model}
    stimuli::Vector{S}
    model::M
    n_rep::Int64=length(stimuli)
    summarizer::Function=mean
    tag::String=""
end

function simulate(p::AvgPattern)
    @showprogress "Estimating excitation pattern" pmap(1:p.n_rep) do idx
        map(p.summarizer, compute(p.model, p.stimuli[idx]))
    end
end

# Define DeltaPattern (pattern emitting summarized difference of two neurograms for pairs of stimuli)
@with_kw struct DeltaPattern{S, M} <: DPattern where {S <: AbstractStimulus, M <: Model}
    stimuli::Vector{Tuple{S, S}}
    model::M
    n_rep::Int64=length(stimuli)
    comparator::Function=-
end

function simulate(p::DeltaPattern)
    @showprogress "Estimating delta pattern" pmap(1:p.n_rep) do idx
        μ₁ = map(p.summarizer, compute(p.model, p.stimuli[idx][1]))
        μ₂ = map(p.summarizer, compute(p.model, p.stimuli[idx][2]))
        p.comparator.(μ₂, μ₁)
    end
end

# Generic function for implementing ID for an EPattern
function id(p::EPattern; accesses=nothing, connector="_", kwargs...)
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
    stim_part = id(p.stimuli[1])

    # Return combination
    return main_part * connector * stim_part
end

# Generic function for implementing ID for a DPattern
function id(p::DPattern; accesses=nothing, connector="_", kwargs...)
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
    stim_part_1 = id(p.stimuli[1][1])
    stim_part_2 = id(p.stimuli[1][2])

    # Return combination
    return main_part * connector * stim_part_1 * stim_part_2
end