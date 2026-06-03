export GammatoneFilterbank

"""
    GammatoneFilterbank <: Model

Auditory gammatone filterbank model.
"""
@with_kw struct GammatoneFilterbank <: Model
    fs::Float64=100e3
    cf::Vector{Float64}=[1000.0]
    coi::Vector{Int64}=1:length(cf)
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)
    bw_factor::Float64=1.0
end

# Some constructors for convenience
GammatoneFilterbank(cf::Float64; kwargs...) = GammatoneFilterbank(; cf=[cf], kwargs...)
GammatoneFilterbank(cf::Vector{Float64}; kwargs...) = GammatoneFilterbank(; cf=cf, kwargs...)

# compute(model, stimulus) maps from stimulus to response
function _compute(m::GammatoneFilterbank, x::Vector{Float64})
    filterbank_gammatone(
        x,
        m.cf,
        m.fs;
        bw_factor=m.bw_factor,
    )
end

function compute(m::GammatoneFilterbank, x::Vector{Float64})
    extract(m, _compute(m, x))
end

islogout(m::GammatoneFilterbank) = true