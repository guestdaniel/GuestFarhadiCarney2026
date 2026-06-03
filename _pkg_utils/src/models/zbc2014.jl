export AuditoryNerveZBC2014

"""
    AuditoryNerveZBC2014 <: Model

Model type for simulating rate from Zilany, Bruce, and Carney (2014) auditory-nerve model
"""
@with_kw struct AuditoryNerveZBC2014 <: Model
    # General parameters
    fs::Float64=100e3
    cf::Vector{Float64}=[1000.0]
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)
    audiogram::Audiogram=Audiogram()

    # IHC parameters
    cohc::Vector{Float64}=fit_audiogram(audiogram, cf)[1]
    cihc::Vector{Float64}=fit_audiogram(audiogram, cf)[2]
    species::String="human"

    # ANF parameters
    fiber_type::String="high"
    power_law::String="approximate"
    fractional::Bool=true
    n_fiber::Int64=1
end

function compute(m::AuditoryNerveZBC2014, x::AbstractVector{Float64}, idx_chan::Int64)
    # Basilar membrane + inner hair cell stage
    ihc_out = sim_ihc_zbc2014(
        x,
        m.cf[idx_chan];
        cohc=m.cohc[idx_chan],
        cihc=m.cihc[idx_chan],
        species=m.species,
    )

    # Auditory nerve stage
    anf_out = zeros(size(ihc_out))
    for _ in 1:m.n_fiber
        anf_out .+= sim_anrate_zbc2014(
            ihc_out,
            m.cf[idx_chan];
            fiber_type=m.fiber_type,
            power_law=m.power_law,
            fractional=m.fractional,
        )
    end
    anf_out = anf_out ./ m.n_fiber
    return anf_out
end

function compute(m::AuditoryNerveZBC2014, x::AbstractVector{Float64})
    map(idx_chan -> compute(m, x, idx_chan), 1:m.n_chan)
end


