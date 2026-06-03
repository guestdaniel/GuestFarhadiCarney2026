export GFC2024_FakeBin

"""
    GFC2024_FakeBin <: Model

(Fake) binaural offshoot of auditory efferent model of Guest, Farhadi, and Carney (2024) 

Variant of the 2024 efferent model with a "fake" binaural system set up to simulate 
contralateral acoustic stimulation experiments. In essence, it consists of two submodels,
an ipsilateral model (`model_ipsi`) and a contralateral model (`model_contra`). To simulated
responses, the contralateral side is first simulated. The contralateral side always has
multiple CFs spanning a sufficient range around the single CF of the ipsilateral side. After
simulating the contralateral response, the gainpostmix signal from the contra channel with 
a CF matching the ipsi CF is used as an input to ipsi model to determine the gain factor.
This, in essence, simulates so-called contralateral acoustic stimulation (CAS) paradigms,
under the key assumption that the ipislateral efferent system does not contribute 
consequentially to responses in the contralateral ear. 

There are a few key notes about the setup that are worth keeping in mind:
    - The CF of the ipsilateral side is always scalar-valued and equal to the middle CF of 
      the contralateral side (selected by `middle(cf)`).
    - The contralateral model has `clip_left=false`, so that its ouptut gain-factor vector
      used in the ipsilateral model has the correct length, taking into account the "extra"
      simulation time before t=0 usually hidden from the user.
    - When the contralateral responses are simulated, the method `compute(::Model, ::Vector)`
      is used. When the ipsilateral responses are simulated, `sim_gfc2023_dict` is directly
      called so that the gain factor can be passed in and `moc_fix_gain=true` can also be
      passed in.
"""
struct GFC2024_FakeBin <: Model
    model_ipsi::GFC2024
    model_contra::GFC2024
end

# Some constructors for convenience
function GFC2024_FakeBin(
    cf::Vector{Float64}=LogRange(0.5e3, 8e3, 31); 
    cf_ipsi=NaN,
    stage="hsr", 
    moc_width=1.5, 
    moc_weight=0.0,
    moc_beta=0.2,
    moc_offset=5.0,
    kwargs...
) 
    # Determine cf_ipsi
    cf_ipsi = isnan(cf_ipsi) ? middle(cf) : cf_ipsi

    # Handle MOC weight (if we get a scalar, make it a vector)
    if !(moc_weight isa Vector)
        moc_weight = fill(moc_weight, length(cf))
    end
    if !(moc_beta isa Vector)
        moc_beta = fill(moc_beta, length(cf))
    end
    if !(moc_offset isa Vector)
        moc_offset = fill(moc_offset, length(cf))
    end
 
    # Create contra model
    model_contra = GFC2024(; 
        cf=cf, 
        stage="gainpostmix", 
        clip_left=false, 
        moc_width=moc_width, 
        moc_weight=moc_weight,
        moc_offset=moc_offset,
        moc_beta=moc_beta,
        coi=[argmin(abs.(cf .- cf_ipsi))],  # make sure gain signal comes from channel closest to cf_ipsi
        kwargs...,
    )

    # Create ipsi model
    model_ipsi = GFC2024(; 
        cf=[isnan(cf_ipsi) ? middle(cf) : cf_ipsi], 
        moc_width=0.0, 
        moc_weight=[0.0],    # prevent any MOC activity
        moc_beta=[1.0],      # arbtirary
        moc_offset=[100.0],  # arbitrary 
        stage=stage, 
        kwargs...,
    )

    # Warn if cf_ipsi is not an exact match to one of the cf values
    if !isnan(cf_ipsi) && !any(isapprox.(cf_ipsi, cf; rtol=0.01))
        cf_match = cf[argmin(abs.(cf .- cf_ipsi))]
        @warn "Ipsilateral CF $(cf_ipsi) Hz does not have exact match in contralateral model, closest match is $(cf_match)!"
    end

    # Combine into single struct
    GFC2024_FakeBin(model_ipsi, model_contra)
end

# Constructor for requesting single CF
function GFC2024_FakeBin(cf::Float64; n_cf=51, species="cat", kwargs...)
    # Branch CF range based on species
    if species == "cat"
        # Cat default goes from 0.2 to 40 (?) kHz  # TODO Look up reasonable values here
        GFC2024_FakeBin(LogRange(0.2e3, 40e3, n_cf); cf_ipsi=cf, species=species, kwargs...)
    else
        # Cat default goes from 0.2 to 20 kHz
        GFC2024_FakeBin(LogRange(0.2e3, 20e3, n_cf); cf_ipsi=cf, species=species, kwargs...)
    end
end

# Some generic methods that need to be extended
samprate(m::GFC2024_FakeBin) = samprate(m.model_ipsi)
ipsi(m::GFC2024_FakeBin) = m.model_ipsi
contra(m::GFC2024_FakeBin) = m.model_contra
getcf(m::GFC2024_FakeBin) = getcf(ipsi(m))[1]
compute(m::GFC2024_FakeBin, x::Vector) = compute(ipsi(m), x)

# compute(model, stimulus) maps from stimulus to response
_compute(m::GFC2024_FakeBin, x::BinauralStimulus) = _compute(m, synthesize(x))

function compute(m::GFC2024_FakeBin, x::Tuple{Vector{Float64}, Vector{Float64}})
    _compute(m, x)[ipsi(m).stage][1]
end

function _compute(m::GFC2024_FakeBin, x::Tuple{Vector{Float64}, Vector{Float64}})
    # Extract stimuli from tuple `x`
    stim_ipsi = x[1]
    stim_contra = x[2]

    # Compute response in contralateral ear
    # Since m.model_contra.stage == "gainpostmix" and m.model_ipsi.coi -> cf, this will be 
    # gain-factor signal for the channel matching the ipsi CF
    g = compute(contra(m), stim_contra)

    # Unpack ipsi model fields
    @unpack_GFC2024 ipsi(m)

    # Compute response in ipsilateral ear with gain fixed
    sim_gfc2023_dict(
        stim_ipsi,
        cf;
        # General parameters
        fs=fs,
        # IHC parameters
        cohc=cohc,
        cihc=cihc,
        species=species,
        # ANF parameters
        fractional=fractional,
        powerlaw_mode=powerlaw_mode,
        # IC parameters
        cn_tau_e=cn_tau_e,
        cn_tau_i=cn_tau_i,
        cn_delay=cn_delay,
        cn_amp=cn_amp,
        cn_inh=cn_inh,
        ic_tau_e=ic_tau_e,
        ic_tau_i=ic_tau_i,
        ic_delay=ic_delay,
        ic_amp=ic_amp,
        ic_inh=ic_inh,
        # MOC params
        moc_cutoff=moc_cutoff,
        moc_beta=moc_beta,
        moc_offset=moc_offset,
        moc_minval=moc_minval,
        moc_maxval=moc_maxval,
        moc_weight=moc_weight,
        moc_width=moc_width,
        moc_delay=moc_delay,
        # Misc params
        dur_pad_left=dur_pad_left,
        clip_left=clip_left,
        dur_pad_right=dur_pad_right,
        clip_right=clip_right,
        # Misc params and gain fixing
        gain=[g],
        moc_fix_gain=true,
    )
end

# custom viz(model, response) function
function viz(m::GFC2024_FakeBin)
    fig = Figure()

    # Plot variation between offset and CF
    ax = Axis(fig[1, 1]; ylabel="Offset (sp/s)", xticklabelsvisible=false, xlabelvisible=false, xscale=log10)
    lines!(ax, getcf(contra(m)) ./ 1e3, contra(m).moc_offset; color=:black)

    # Plot variation between β and CF
    ax = Axis(fig[2, 1]; xlabel="CF (kHz)", ylabel="β", xscale=log10)
    lines!(ax, getcf(contra(m)) ./ 1e3, contra(m).moc_beta; color=:black)

    fig
end

viz(m::GFC2024_FakeBin, s::AbstractStimulus; kwargs...) = viz(m, synthesize(s); kwargs...)

function viz(m::GFC2024_FakeBin, x::Tuple{Vector{Float64}, Vector{Float64}}; kwargs...)
    viz(m, x, _compute(m, x); kwargs...)
end

function viz(
    model::GFC2024_FakeBin, 
    x::Tuple{Vector{Float64}, Vector{Float64}},
    resp::Dict;
    stages=["ihc", "hsr", "lsr", "gainpostmix"],
)
    # Unsplat stimulus `x`
    stim_ipsi = x[1]
    stim_contra = x[2]

    # Create figure and plot each element
    fig = Figure(; size=(600, max(500, 100*length(stages))))
    axs = [Axis(fig[i, 1]) for i in 1:(length(stages)+2)]
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    lines!(axs[1], timevec(stim_ipsi, samprate(model)), stim_ipsi)
    lines!(axs[2], timevec(stim_contra, samprate(model)), stim_contra)
    for (idx, key) in enumerate(stages)
        r = resp[key][1]
        lines!(axs[idx+2], timevec(r, samprate(model)), r)
    end
    xlims!.(axs, (extrema(timevec(stim_ipsi, samprate(model))) .+ (0.0, 0.025))...)
    for (idx, name) in enumerate(vcat("Ipsi stimulus", "Contra stimulus", stages))
        Label(fig[idx, 2], name; tellheight=false)
    end
    fig
end

# custom neurogram(model, response) function 
# this custom function shows the contralateral neurogram only!
function neurogram(model::GFC2024_FakeBin, stim::BinauralStimulus)
    neurogram(contra(model), contra(stim))
end