export NHPPSpikingNeuron

struct NHPPSpikingNeuron{M} <: Model where {M<:Model}
    model::M
    n_rep::Int64
    binsize::Float64
    pretransform::Function
end

NHPPSpikingNeuron(model; n_rep=100, binsize=1e-3, pretransform=identity) = NHPPSpikingNeuron(model, n_rep, binsize, pretransform)

getcoi(model::NHPPSpikingNeuron) = getcoi(model.model)

# Method to compute single response for single or multichannel models
function compute(model::NHPPSpikingNeuron{M}, stim::AbstractStimulus) where {M<:Model}
    # Compute response and pretransform
    R = compute(model.model, stim)  # Vector{Float64} for singlechannel, Vector{Vector{Float64}} for multichannel
    R = model.pretransform(R)

    # Figure out what T is (assuming we always left clip)
    T = model.model.clip_right ? dur(stim) : dur(stim) + model.model.dur_pad_right

    # Branch based on whether this is a multichannel or singlechannel model
    if issinglechannel(model)
        # For singlechannel model transform Vector{Float64} R into Vector{Vector{Float64}} 
        # of n_rep spiketimes by sampling from NHPP
        S = sample(model, R)
        x = calc_psth(S, T, model.binsize)

        # Interpolate PSTH to match resolution of simulation
        t = 0.0:(1/samprate(stim)):(T-1/samprate(stim))
        binpos = 0.5 .* (x.edges[1][2:end] .+ x.edges[1][1:(end-1)])
        itp = Interpolations.interpolate((binpos,), Float64.(x.weights), Gridded(Constant()))
        etp = extrapolate(itp, Interpolations.Flat())
        etp.(t)
    else
        # For multichannel model transform Vector{Vector{Float64}} R into 
        # Vector{Vector{Vector{Float64}}} of n_rep spiketimes by sampling from NHPP
        map(R) do chan
            S = sample(model, chan)
            x = calc_psth(S, T, model.binsize)
            t = 0.0:(1/samprate(stim)):(T-1/samprate(stim))
            binpos = 0.5 .* (x.edges[1][2:end] .+ x.edges[1][1:(end-1)])
            itp = Interpolations.interpolate((binpos,), Float64.(x.weights), Gridded(Constant()))
            etp = extrapolate(itp, Interpolations.Flat())
            etp.(t)
        end
    end
end

# Method to compute multiple responses in-place for single or multichannel models
function compute(model::NHPPSpikingNeuron{M}, stim::AbstractStimulus, N::Int64) where {M<:Model}
    # Compute response and pretransform
    Rs = compute(model.model, stim, N)  # Vector{Float64} for singlechannel, Vector{Vector{Float64}} for multichannel

    # Map over each R
    map(Rs) do R
        R = model.pretransform(R)

        # Figure out what T is (assuming we always left clip)
        T = model.model.clip_right ? dur(stim) : dur(stim) + model.model.dur_pad_right

        # Branch based on whether this is a multichannel or singlechannel model
        if issinglechannel(model)
            # For singlechannel model transform Vector{Float64} R into Vector{Vector{Float64}} 
            # of n_rep spiketimes by sampling from NHPP
            S = sample(model, R)
            x = calc_psth(S, T, model.binsize)

            # Interpolate PSTH to match resolution of simulation
            t = 0.0:(1/samprate(stim)):(T-1/samprate(stim))
            binpos = 0.5 .* (x.edges[1][2:end] .+ x.edges[1][1:(end-1)])
            itp = Interpolations.interpolate((binpos,), Float64.(x.weights), Gridded(Constant()))
            etp = extrapolate(itp, Interpolations.Flat())
            etp.(t)
        else
            # For multichannel model transform Vector{Vector{Float64}} R into 
            # Vector{Vector{Vector{Float64}}} of n_rep spiketimes by sampling from NHPP
            map(R) do chan
                S = sample(model, chan)
                x = calc_psth(S, T, model.binsize)
                t = 0.0:(1/samprate(stim)):(T-1/samprate(stim))
                binpos = 0.5 .* (x.edges[1][2:end] .+ x.edges[1][1:(end-1)])
                itp = Interpolations.interpolate((binpos,), Float64.(x.weights), Gridded(Constant()))
                etp = extrapolate(itp, Interpolations.Flat())
                etp.(t)
            end
        end
    end
end

function sample(model::NHPPSpikingNeuron, λ::Vector{Float64})
    map(x -> nhpp_thinning(λ, samprate(model)), 1:model.n_rep)
end

function sample(model::NHPPSpikingNeuron, λ::Vector{Vector{Float64}})
    map(eachλ -> map(x -> nhpp_thinning(eachλ, samprate(model)), 1:model.n_rep), λ)
end

function viz(
    model::NHPPSpikingNeuron,
    x::Vector{<:AbstractStimulus};
    colors=colorschemes[:Dark2_8],
    linewidth=0.5,
)
    # Compute response
    resps = map(x) do _x
        compute(model, _x)
    end

    # Create figure and plot each element
    fig = Figure(; size=(600, 400))
    axs = [Axis(fig[i, 1]) for i in 1:2]
    hidexdecorations!.(axs[1:(end-1)], ticks=false, grid=false)
    map(zip(x, colors)) do (_x, c)
        w = synthesize(_x)
        lines!(axs[1], timevec(w, samprate(model)), w; color=c)
    end
    map(zip(resps, colors)) do (resp, c)
        lines!(axs[2], timevec(resp, samprate(model)), resp; color=c, linewidth=linewidth)
    end
    xlims!.(axs, (extrema(vcat(timevec.(x)...)) .+ (0.0, 0.025))...)
    ylims!(axs[2], standardylims(getstage(model))...)
    fig
end

viz(model::NHPPSpikingNeuron, stim::AbstractStimulus) = viz(model, [stim])

samprate(model::NHPPSpikingNeuron) = samprate(model.model)
getstage(model::NHPPSpikingNeuron) = model.model.stage
getcf(model::NHPPSpikingNeuron) = getcf(model.model)

