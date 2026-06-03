module Utilities

using AuditorySignalUtils
using CairoMakie
using Colors
using ColorSchemes
using Dates
using DrWatson
using Distributed
using Interpolations
using MAT
using SHA
using Parameters
using ProgressMeter
using DSP
using FFTW
using Crayons
using Statistics
using OnlineStats
using Distributions
using Suppressor
using Random
using Match
using Helios
using CubicSplines
using ZilanyBruceCarney2014
using StatsBase
using SkipNan
using Optim
using Loess
using Interpolations
using Printf
using Format

# Misc/utility functions
include(joinpath("misc", "misc.jl"))
include(joinpath("misc", "smoothing.jl"))

# Core code
include(joinpath("core", "utils.jl"))
include(joinpath("core", "components.jl"))
include(joinpath("core", "audiograms.jl"))
include(joinpath("simulations", "configs.jl"))
include(joinpath("core", "memo.jl"))
include(joinpath("core", "plots.jl"))

# Stimulus code
include(joinpath("stimuli", "filtering.jl"))
include(joinpath("stimuli", "stimuli.jl"))
include(joinpath("stimuli", "tones.jl"))
include(joinpath("stimuli", "noises.jl"))
include(joinpath("stimuli", "compound.jl"))
include(joinpath("stimuli", "clicks.jl"))

# Modeling code
include(joinpath("models", "models.jl"))
include(joinpath("models", "zbc2014.jl"))
include(joinpath("models", "nc2004.jl"))
include(joinpath("models", "gfc2023.jl"))
include(joinpath("models", "gfc2024.jl"))
include(joinpath("models", "gfc2024_fakebin.jl"))
include(joinpath("models", "spiking.jl"))
include(joinpath("models", "gammatone.jl"))

# Simulations and Experiments
include(joinpath("simulations", "simulations.jl"))
include(joinpath("simulations", "rlf.jl"))
include(joinpath("simulations", "tuning_curves.jl"))
include(joinpath("simulations", "misc.jl"))
include(joinpath("simulations", "crlb.jl"))
include(joinpath("simulations", "crlb_abstraction_X.jl"))

# Set up standards for types
include(joinpath("standards", "standards.jl"))

end  # module
