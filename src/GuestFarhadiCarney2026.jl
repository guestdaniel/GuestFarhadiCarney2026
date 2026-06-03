module GuestFarhadiCarney2026

using ADTypes
using AuditorySignalUtils
using CairoMakie
using CSV
using Colors
using ColorSchemes
using DataFrames
using DataFramesMeta
using Distributed
using Distributions
using DrWatson
using DSP
using FFTW
using Helios
using Interpolations
using Loess
using Match
using ProgressMeter
using Parameters
using SkipNan
using StatsBase
using Statistics
using Utilities
using ZilanyBruceCarney2014

# Include
include("PrecursorRLF.jl")                                   
include(joinpath("WarrenLiberman", "WarrenLiberman1989.jl"))  # Code for Warren CAS
include("utils.jl")
include("stimuli.jl")
include("simulations.jl")
include("figures.jl")
include("supp_figures.jl")

end