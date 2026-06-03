export GFC2023

"""
    GFC2023 <: Model

In-progress auditory model, Guest, Farhadi, and Carney (2023) 
"""
@with_kw struct GFC2023 <: Model
    # General parameters
    fs::Float64=100e3
    cf::Vector{Float64}=[1000.0]
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)
    audiogram::Audiogram=Audiogram()
    stage::String="hsr"

    # IHC parameters
    cohc::Vector{Float64}=fit_audiogram(audiogram, cf)[1]
    cihc::Vector{Float64}=fit_audiogram(audiogram, cf)[2]
    species::String="human"

    # ANF parameters
    powerlaw_mode::Int64=2
    fractional::Bool=false

    # MOC parameters
    moc_weight_ic::Float64=8.0
    moc_weight_wdr::Float64=2.0
    moc_width_wdr::Float64=0.5
    moc_offset_wdr::Float64=0.0
    moc_offset_ic::Float64=0.0
end

GFC2023(cf::Float64; kwargs...) = GFC2023(; cf=[cf], kwargs...)
GFC2023(cf::Vector{Float64}; kwargs...) = GFC2023(; cf=cf, kwargs...)

function compute(m::GFC2023, x::Vector{Float64})
    sim_gfc2023_dict(
        x,
        m.cf;
        cohc=m.cohc,
        cihc=m.cihc,
        species=m.species,
        fractional=m.fractional,
        powerlaw_mode=m.powerlaw_mode,
        moc_weight_ic=m.moc_weight_ic,
        moc_weight_wdr=m.moc_weight_wdr,
        moc_width_wdr=m.moc_width_wdr,
        moc_offset_wdr=m.moc_offset_wdr,
        moc_offset_ic=m.moc_offset_ic,
        dur_pad_left=0.02,
        clip_left=true,
        dur_pad_right=0.015,
        clip_right=false,
    )[m.stage]
end

function compute(m::GFC2023, x::Vector{Float64}, stages)
    resp = sim_gfc2023_dict(
        x,
        m.cf;
        cohc=m.cohc,
        cihc=m.cihc,
        species=m.species,
        fractional=m.fractional,
        powerlaw_mode=m.powerlaw_mode,
        moc_weight_ic=m.moc_weight_ic,
        moc_weight_wdr=m.moc_weight_wdr,
        moc_width_wdr=m.moc_width_wdr,
        moc_offset_wdr=m.moc_offset_wdr,
        moc_offset_ic=m.moc_offset_ic,
        dur_pad_left=0.01,
        clip_left=true,
        dur_pad_right=0.015,
        clip_right=false,
    )
    [resp[stage] for stage in stages]
end
