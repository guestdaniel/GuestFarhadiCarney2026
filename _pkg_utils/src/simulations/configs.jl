export Config, Default, setrng

"""
    Config

Encodes metadata about simulations and governs their metabehavior (e.g., caching)

Configs encode a variety of configuration settings (e.g., cache settings), config (e.g.,
username and machine), and other meta-information that is useful in controlling how
simulations are carried out, reconstructing simulations in the future, and minimizing
unnecessary recomputation.

The current available Configs are:
- `Config`: Standard Config for most daily driving, governs caching behavior and stores
  information about time, Git status, and machines.
"""
abstract type Config end

@with_kw struct Default <: Config 
    # Misc info
    codename::String="helloworld"
    # Git info
    gitcommit::Union{Nothing, String}=@suppress DrWatson.gitdescribe(projectdir())
    # Machine info
    machine::String=(Sys.KERNEL == :Linux) ? ENV["HOSTNAME"] : ENV["COMPUTERNAME"]
    user::String=(Sys.KERNEL == :Linux) ? ENV["USER"] : ENV["USERNAME"]
    platform::String=Sys.MACHINE
    # Seed info
    rng::Xoshiro=copy(Random.default_rng())
    seed::Int64=(round(rand() * 100_000_00))
    # Paths
    path_cache::String=(Sys.KERNEL == :Linux) ? "/scratch/dguest2/cl_cache" : "\\home\\daniel\\cl_cache"
    path_out::String=(Sys.KERNEL == :Linux) ? "/scratch/dguest2/cl_fig" : "\\home\\daniel\\cl_fig"
    # Cache resolution flags
    resolve_codename::Bool=false
    resolve_commit::Bool=false
    resolve_machine::Bool=false
    resolve_rng::Bool=false
    # Cache behavior flags
    load_from_cache::Bool=true
    save_to_cache::Bool=true
    # Display flags
    disp_progress::Bool=true
end

Config(args...; kwargs...) = Default(args...; kwargs...)

function setrng(config::Config)
    Random.seed!(config.seed)
end

function Base.display(config::Config)
    print(
        """
        $(string(typeof(config))) config with codename: $(config.codename)
        Running on $(config.machine) / $(config.platform) under user $(config.user) 
        Name resolution:
            codename: $(config.resolve_codename) | machine: $(config.resolve_machine) | commit: $(config.resolve_commit) | rng: $(config.resolve_rng)
        Caching behavior:
            load: $(config.load_from_cache) | save: $(config.save_to_cache)
        """
    )
end

cachepath(config::Config) = config.path_cache

function id(config::Config; connector="_")
    # Figure out what information to encode in id based on flags
    accesses = Symbol[]
    # If resolve_commit, we include information about the git commit
    if config.resolve_commit
        push!(accesses, :gitcommit)
    end
    # If resolve_machine, we include information about the machine, user, and platform
    if config.resolve_machine
        push!(accesses, :machine)
        push!(accesses, :user)
        push!(accesses, :platform)
    end
    # If resolve_rng, we include the Seed
    if config.resolve_rng
        push!(accesses, :rng)
        push!(accesses, :seed)
    end
    # If resolve_codename, we include the codename
    if config.resolve_codename
        push!(accesses, :codename)
    end

    # Construct id using DrWatson.savename
    savename(
        config;
        accesses=accesses,
        allowedtypes=(Real, String, Symbol, DateTime, Function, Xoshiro),
        connector=connector,
    )
end

Base.hash(config::Config) = Base.hash(id(config))