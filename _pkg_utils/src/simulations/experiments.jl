export Experiment           # abstract tpyes
#export Subexperiments
export status, list, setup, purge, viz       # functions

"""
    Experiment

Abstract type for bundling together related code for sequences or groups of simulations
"""
abstract type Experiment end

# Declare methods that must be provided by user
function setup(::Experiment) end
function viz(::Experiment; kwargs...) end
function list(::Experiment) end
function purge(::Experiment) end

# Declare methods provided by us
function Base.run(experiment::Experiment; config::Config=Default())
    # Set up simulations and context
    sims = setup(experiment)

    # Loop through simulations and evaluate
    @info "Running $experiment"
    @showprogress for sim in sims
        @memo config simulate(sim)
    end

    # Visualize
    viz(experiment; config=config)
end

function status(experiment::Experiment, func::Function=simulate; config::Config=Default())
    @info "Checking status of $experiment"
    display(config)
    sims = setup(experiment)
    cached = map(sims) do sim
        isfile(cachepath(config, func, sim))
    end
    if allcached(cached)
        @info "Experiment complete and cached!"
    else
        @info "Experiment incomplete!"
        display(cached)
    end
    return cached
end

# function purge(experiment::Experiment, config::Config=Default())
#     @info "Purging cache for $experiment"
#     sims = setup(experiment)
#     map(sims) do sim
#         rm(cachepath(config, sim))
#     end
#     @info "Cache purged!"
# end

allcached(cached::Any) = all(map(allcached, cached))
allcached(cached::Vector{Bool}) = all(cached)
allcached(cached::Matrix{Bool}) = all(cached)
sumcache(cached::Any) = sum(map(sumcache, cached))
sumcache(cached::Vector{Bool}) = sum(cached)
sumcache(cached::Matrix{Bool}) = sum(cached)
lencache(cached::Any) = sum(map(lencache, cached))
lencache(cached::Vector{Bool}) = length(cached)
lencache(cached::Matrix{Bool}) = length(cached)