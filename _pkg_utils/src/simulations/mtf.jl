export NoiseMTF, ToneMTF, axis

"""
    NoiseMTF

Type for generating responses to vector of SAM noises for given model
"""
@with_kw struct NoiseMTF{M} <: Simulation where {M <: Model}
    stimuli::Vector{SAMNoise}
    model::M
    summary::Function=mean
    fms::Vector{Float64}=[stim.fm for stim in stimuli]
    n_fm::Float64=length(fms)
    fm_low::Float64=minimum(fms)
    fm_high::Float64=maximum(fms)
    level::Float64=stimuli[1].l
end

function NoiseMTF(model::M, fm_low=8.0, fm_high=512.0, n_fm=14; fs=100e3, kwargs...) where {M <: Model}
    stimuli = vcat(
        SAMNoise(; fm=0.0, fs=fs, kwargs...),
        map(fm -> SAMNoise(; fm=fm, fs=fs, kwargs...), LogRange(fm_low, fm_high, n_fm)),
    )
    NoiseMTF(; stimuli=stimuli, model=model)
end

function simulate(s::NoiseMTF) 
    map(s.stimuli) do stim
        s.summary(extract(s.model, compute(s.model, stim)))
    end
end

function axis(s::NoiseMTF)
    [stim.fm for stim in s.stimuli]
end