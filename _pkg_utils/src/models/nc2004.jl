export InferiorColliculusSFIE, InferiorColliculusSFIEBE, InferiorColliculusSFIEBS, InferiorColliculusSFIE_Multiunit, sfie

# First, declare a few functions that we'll need to actually compute the inferior-colliculus
# model responses, including:
# - `get_α_normalized`: returns coefficients for an SFIE "alpha" filter
# - `sim_sfie_nc2004`: returns output given input and single-stage set of SFIE parameters

"""
    get_α_normalized(τ, fs, t)

Get coefficients for a normalized "alpha" filter in SFIE model
"""
function get_α_normalized(τ, fs, t)
    α = exp(-1.0 / (fs * τ))
    norm = 1 ./ (τ^2 .* (exp(-t/τ) .* (-t/τ-1) + 1))
    B = [0, α]
    A = [1, -2*α, α^2] * fs * 1 ./ norm
    return B, A
end

"""
    sim_sfie_nc2004(input)
    sim_sfie_nc2004(input[; τ_e=0.5e-3, τ_i=2.0e-3, d_i=1e-3, S=0.6, A=1.5, fs=100e3])

Given input and SFIE parameter set, simulate output rate
"""
function sim_sfie_nc2004(
    input::AbstractVector{Float64};
    τ_e=0.5e-3,
    τ_i=2.0e-3,
    d_i=1.0e-3,
    S=0.6,
    A=1.5,
    fs=100e3,
)
    # Delay
    E = input
    I = shiftsignal(input, Int(floor(d_i * fs)))

    # Filter with frequency-domain α functions
    (b, a) = get_α_normalized(τ_e, fs, 1)
    E = filt(b, a, E)
    @. E = E * A * (1.0/fs)
    (b, a) = get_α_normalized(τ_i, fs, 1)
    I = filt(b, a, I)
    @. I = I * A * S * (1.0/fs)

    # Return output
    output = ((E .- I) .+ abs.(E .- I))./2
    return output
end

# Next, we declare abstract type for all SFIE models (InferiorColliculusSFIE), as well as
# concrete subtypes for BE and BS cell models. We also provide a method called `sfie` that
# has a separate method for BE and BS cell types. The obligatory `compute` method for
# subtypes of `Model` is provided at the level of `InferiorColliculusSFIE` and calls `sfie`. 

"""
    InferiorColliculusSFIE <: Model 

Abstract type for SFIE model units
"""
abstract type InferiorColliculusSFIE <: Model end

function compute(m::InferiorColliculusSFIE, x::AbstractVector{Float64})
    map(chan -> sfie(m, chan), compute(m.frontend, x))
end

"""
    InferiorColliculusSFIEBE <: InferiorColliculusSFIE

Model type for band-enhanced inferior-colliculus model neurons
"""
@with_kw struct InferiorColliculusSFIEBE <: InferiorColliculusSFIE
    # Model interface
    frontend::AuditoryNerveZBC2014=AuditoryNerveZBC2014()
    fs::Float64=frontend.fs
    cf::Vector{Float64}=frontend.cf
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)

    # Cochlear nucleus parameters
    τ_e_cn::Float64=0.5e-3
    τ_i_cn::Float64=2e-3
    d_i_cn::Float64=1e-3
    S_cn::Float64=0.6
    A_cn::Float64=1.5

    # Inferior colliculus parameters
    τ_e_ic::Float64=1e-3
    τ_i_ic::Float64=1.5e-3
    d_i_ic::Float64=2e-3
    S_ic::Float64=0.9
    A_ic::Float64=1.0
end

function sfie(m::InferiorColliculusSFIEBE, x::AbstractVector{Float64})
    # Cochlear nucleus stage
    cn_out = sim_sfie_nc2004(
        x;
        τ_e=m.τ_e_cn,
        τ_i=m.τ_i_cn,
        d_i=m.d_i_cn,
        S=m.S_cn,
        A=m.A_cn,
    )
    # Inferior colliculus stage
    ic_out = sim_sfie_nc2004(
        cn_out;
        τ_e=m.τ_e_ic,
        τ_i=m.τ_i_ic,
        d_i=m.d_i_ic,
        S=m.S_ic,
        A=m.A_ic,
    )
    return ic_out
end

"""
    InferiorColliculusSFIEBS <: InferiorColliculusSFIE

Model type for band-suppressed inferior-colliculus model neurons
"""
@with_kw struct InferiorColliculusSFIEBS <: InferiorColliculusSFIE
    # Model interface
    frontend::AuditoryNerveZBC2014=AuditoryNerveZBC2014()
    fs::Float64=frontend.fs
    cf::Vector{Float64}=frontend.cf
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)

    # Cochlear nucleus parameters
    τ_e_cn::Float64=0.5e-3
    τ_i_cn::Float64=2e-3
    d_i_cn::Float64=1e-3
    S_cn::Float64=0.6
    A_cn::Float64=1.5

    # Inferior colliculus (BE) parameters
    τ_e_ic::Float64=1e-3
    τ_i_ic::Float64=1.5e-3
    d_i_ic::Float64=2e-3
    S_ic::Float64=0.9
    A_ic::Float64=1.0

    # Inferior colliculus (BS) parameters
    τ_i_bs::Float64=τ_i_ic
    d_i_bs::Float64=1.0e-3
    S_bs::Float64=4.0
    A_bs::Float64=0.2
end

function sfie(m::InferiorColliculusSFIEBS, x::AbstractVector{Float64})
    # Cochlear nucleus stage
    cn_out = sim_sfie_nc2004(
        x;
        τ_e=m.τ_e_cn,
        τ_i=m.τ_i_cn,
        d_i=m.d_i_cn,
        S=m.S_cn,
        A=m.A_cn,
    )

    # BE stage
    E = cn_out
    I = shiftsignal(cn_out, Int(floor(m.d_i_ic * m.fs)))
    (b, a) = get_α_normalized(m.τ_e_ic, m.fs, 1)
    E = filt(b, a, E)
    @. E = E * m.A_ic * (1.0/m.fs)
    (b, a) = get_α_normalized(m.τ_i_ic, m.fs, 1)
    I = filt(b, a, I)
    @. I = I * m.A_ic * m.S_ic * (1.0/m.fs)
    be_out = ((E .- I) .+ abs.(E .- I))./2

    # BS stage
    I = shiftsignal(be_out, Int(floor(m.d_i_bs * m.fs)))
    (b, a) = get_α_normalized(m.τ_i_bs, m.fs, 1)
    I = filt(b, a, I)
    @. I = I * m.A_bs * m.S_bs * (1.0/m.fs)
    bs_out = ((E .- I) .+ abs.(E .- I))./2

    return bs_out
end

"""
    InferiorColliculusSFIE_Multiunit <: InferiorColliculusSFIE

Model type for multiple BE/BS units combined together
"""
@with_kw struct InferiorColliculusSFIE_Multiunit <: InferiorColliculusSFIE
    # Model interface
    frontend::AuditoryNerveZBC2014=AuditoryNerveZBC2014()
    fs::Float64=frontend.fs
    cf::Vector{Float64}=frontend.cf
    n_chan::Int64=length(cf)
    cf_low::Float64=minimum(cf)
    cf_high::Float64=maximum(cf)

    # Cochlear nucleus parameters
    τ_e_cn::Float64=0.5e-3
    τ_i_cn::Float64=2e-3
    d_i_cn::Float64=1e-3
    S_cn::Float64=0.6
    A_cn::Float64=1.5

    # IC parameters (represented as multiple BE and multiple BS stages)
    units_be::Vector{InferiorColliculusSFIEBE}=InferiorColliculusSFIEBE[]
    units_bs::Vector{InferiorColliculusSFIEBS}=InferiorColliculusSFIEBS[]
end

function InferiorColliculusSFIE_Multiunit(
    data; 
    kwargs...
)
    # Each element in data is a tuple of a model DataType and parameters, we need to turn
    # them into models
    units_be = InferiorColliculusSFIEBE[]
    units_bs = InferiorColliculusSFIEBS[]
    for (model, param) in data
        if model == InferiorColliculusSFIEBE
            push!(units_be, model(; kwargs..., param...))
        else
            push!(units_bs, model(; kwargs..., param...))
        end
    end

    # Create InferiorColliculusSFIE_Multiunit
    InferiorColliculusSFIE_Multiunit(; units_be=units_be, units_bs=units_bs, kwargs...)
end

function compute(m::InferiorColliculusSFIE_Multiunit, x::AbstractVector{Float64})
    # Collect frontend responses
    resps_frontend = compute(m.frontend, x)

    # Loop over each BE/BS unit in the set and simulate, then concatenate
    be_rates = map(m.units_be) do unit
        map(chan -> sfie(unit, chan), resps_frontend)
    end

    bs_rates = map(m.units_bs) do unit
        map(chan -> sfie(unit, chan), resps_frontend)
    end

    vcat(vcat(be_rates)..., vcat(bs_rates)...)
end

function id(comp::InferiorColliculusSFIE_Multiunit; accesses=nothing, connector="_", kwargs...)
    main_part = savename(
        string(typeof(comp)),
        comp; 
        accesses=accesses === nothing ? fieldnames(typeof(comp)) : accesses,
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

    be_part = join(map(id, comp.units_be), connector)
    bs_part = join(map(id, comp.units_bs), connector)

    join([main_part, be_part, bs_part], connector)
end

