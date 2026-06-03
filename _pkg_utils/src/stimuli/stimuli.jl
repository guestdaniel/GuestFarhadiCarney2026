export AbstractStimulus, AbstractNoiseStimulus, AbstractCompoundStimulus, AbstractSINStimulus, Silence
export synthesize, parts, freq, level, dur, onset, onsets, offset, offsets, timeaxis, samprate, sampint, viz, idxswin, bw, freqrange, snr, gettitle, modfreq, vector_strength, syncrate, Δrate

# Define type hierarchy and methods (excluding visualization methods) for stimuli
#
# All stimuli are subtypes of AbstractStimulus
# Top-level methods are defined for this type, relying on presence of `level`, `dur`, and
# `fs` attributes. Methods `onset`, `offset`, `onsets`, `offsets`, an `idxswin` may need to 
# overwritten by concrete types.
abstract type AbstractStimulus <: Component end
function synthesize(::AbstractStimulus) end
AuditorySignalUtils.samples(x::AbstractStimulus) = samples(dur(x), samprate(x)) 
timeaxis(x::AbstractStimulus) = 0.0:sampint(x):nextfloat(dur(x) - sampint(x))
freqaxis(x::AbstractStimulus, N=1) = LinRange(0.0, samprate(x), samples(x)*N)
level(x::AbstractStimulus) = x.level
dur(x::AbstractStimulus) = x.dur
samprate(x::AbstractStimulus) = x.fs
sampint(x::AbstractStimulus) = 1/samprate(x)
onset(x::AbstractStimulus) = 0.0
offset(x::AbstractStimulus) = dur(x) - sampint(x)
onsets(x::AbstractStimulus) = [onset(x)]
offsets(x::AbstractStimulus) = [offset(x)]
idxswin(x::AbstractStimulus) = sampleat(onset(x), samprate(x)):sampleat(offset(x), samprate(x))
gettitle(x::AbstractStimulus) = string(typeof(x))
modfreq(x::AbstractStimulus) = NaN
Base.length(x::AbstractStimulus) = samples(dur(x), samprate(x))
AuditorySignalUtils.timevec(x::AbstractStimulus) = timevec(dur(x), samprate(x))

# Define methods for common analysis stimuli (`mean`, `rms`, `maximum`, `minimum`) that
# include a second argument of an AbstractStmulus subtype. This is a convenient way of 
# creating analysis functions that can incorporate information about the stimulus parameters
# into the analysis. We define below some versions for a Vector containing float data 
# (implicitly, single channel response vector) and also versions for Vectors containing
# other Vectors (implicitly, multichannel response vector).
Statistics.mean(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = mean(r, args...)
Statistics.std(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = std(r, args...)
Base.maximum(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = maximum(r, args...)
Base.minimum(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = minimum(r, args...)
DSP.rms(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = rms(r, args...)
vector_strength(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = synccoef_cont(r, modfreq(x), samprate(x))
Δrate(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = maximum(r) - minimum(r)
syncrate(r::Vector{<:AbstractFloat}, x::AbstractStimulus, args...) = vector_strength(r, x, args...)*mean(r, x, args...)

# Tonal stimuli with a well-defined frequency or F0 are subtypes of AbstractToneStimulus
# Tone-specific methods, such as `freq`, are defined for this type. These rely currently on 
# the existence of a `freq` attribute.
abstract type AbstractToneStimulus <: AbstractStimulus end
freq(x::AbstractToneStimulus) = x.freq


# Noise stimuli without a well-defined frequency or F0 are subtypes of AbstractNoiseStimulus
# Noise-specific methods, such as `sl` and `bw` are defined for this type. These currently
# rely on the existence of `freq_high` and `freq_low` attributes.
abstract type AbstractNoiseStimulus <: AbstractStimulus end
bw(x::AbstractNoiseStimulus) = x.freq_high - x.freq_low
freqrange(x::AbstractNoiseStimulus) = (x.freq_low, x.freq_high)


# Signal-in-noise stimuli are common enough that it is useful to have a specific abstract
# type for them. These methods rely on `signal` and `noise` attributes. At the abstract
# level, we assume `signal` is a concrete subtype of AbstractToneStimulus and `noise` is a
# concrete subtype of AbstractNoiseStimulus. By default, we assume that when we want to mix
# the two synthesized waveforms together that an `idxswin` method is defined for the type 
# and the window lenght should match the length of the signal; the signal is mixed into the
# noise/masker at those samples. Moreover, we assume generally that the signal has a shorter
# duration than the noise, therefore `dur` points to noise duration, and that both signals
# have the same sampling rtae, therefore `samprate` points to the signal sampling rate.
abstract type AbstractSINStimulus <: AbstractStimulus end
freq(x::AbstractSINStimulus) = freq(x.signal)
bw(x::AbstractSINStimulus) = bw(x.noise)
freqrange(x::AbstractSINStimulus) = freqrange(x.noise)
dur(x::AbstractSINStimulus) = dur(x.noise)
samprate(x::AbstractSINStimulus) = samprate(x.signal)
snr(x::AbstractSINStimulus) = level(x.signal) - level(x.noise)

function synthesize(stim::AbstractSINStimulus)
    # Synthesize signal and noise
    signal = synthesize(stim.signal)
    noise = synthesize(stim.noise)

    # Mix and return
    noise[idxswin(stim)] .+= signal
    noise 
end


# A final abstract type is provided for combinations of other stimulus types; this is called
# AbstractCompoundStimulus
abstract type AbstractCompoundStimulus <: AbstractStimulus end
function parts(::AbstractCompoundStimulus) end

# Define a few concrete stimulus types
"""
    Silence <: AbstractCompoundStimulus

Silence with duration `d` and sample rate `fs`
"""
@with_kw struct Silence <: AbstractCompoundStimulus
    dur::Float64=1.0
    fs::Float64=100e3
end

synthesize(stim::Silence) = zeros(samples(stim.dur, stim.fs))


# Visualization methods for stimuli; plot time-domain waveforms and spectra. We first define
# a standard method that works on any stimulus, plotting the time-domain waveform and spectrum
# of the output of `synthesize` side by side.
function viz(
    x::AbstractStimulus; 
    xscale=log10, 
    ylims=(-50, 10), 
    xlims=(0.01, 20), 
    size=(800, 300), 
    highlight=true, 
    fig = Figure(; size=size),
    ax_w = Axis(fig[1, 1]),
    ax_s = Axis(fig[1, 2]; xscale=xscale),
    norm=NaN,
    show_title=false,
    times=[],
    kwargs...
)
    # Optionally add highlight
    if highlight
        if isempty(times)
            vspan!(ax_w, timeat.(extrema(idxswin(x)), samprate(x))...; color=(:red, 0.1))
        else
            vspan!(ax_w, times...; color=(:red, 0.1))
        end
    end

    # Create time waveform and plot
    w = synthesize(x)
    lines!(ax_w, timeaxis(x), w)
    ylims!(ax_w, (1.5 .* extrema(w))...)
    if onsets(x) isa StepRangeLen
        vlines!(ax_w, onsets(x); color=:black, linestyle=:dash)
        vlines!(ax_w, offsets(x); color=:black, linestyle=:dot)
    end

    # Calculate spectrum and plot
    S = 20 .* log10.(abs.(fft(w)))
    lines!(ax_s, freqaxis(x) ./ 1000, S .- (isnan(norm) ? maximum(S) : norm))
    xlims!(ax_s, xlims...)
    ylims!(ax_s, ylims...)

    # Add labels
    ax_w.xlabel = "Time (s)"
    ax_w.ylabel = "Amplitude (Pa)"
    ax_s.xlabel = "Frequency (kHz)"
    ax_s.ylabel = "Level (dB re: max)"

    # Add title
    if show_title
        Label(fig[0, :], gettitle(x); fontsize=20.0)
    end

    fig
end

# Here we define another `viz` method but specifically for AbstractSINStimulus
function viz(
    x::AbstractSINStimulus; 
    xscale=log10, 
    ylims=(-50, 10), 
    xlims=(0.01, 20), 
    size=(900, 300), 
    highlight=true, 
    fig = Figure(; size=size),
    ax_w = Axis(fig[1, 1]),
    ax_s = Axis(fig[1, 2]; xscale=xscale),
    ax_s2 = Axis(fig[1, 3]; xscale=xscale),
    norm=NaN,
    kwargs...
)
    # Optionally add highlight
    if highlight
        vspan!(ax_w, timeat.(extrema(idxswin(x)), samprate(x))...; color=(:red, 0.1))
    end

    # Create time waveform and plot
    w = synthesize(x)
    lines!(ax_w, timeaxis(x), w)
    if onsets(x) isa StepRangeLen
        vlines!(ax_w, onsets(x); color=:black, linestyle=:dash)
        vlines!(ax_w, offsets(x); color=:black, linestyle=:dot)
    end

    # Calculate spectrum and plot
    w_noi = synthesize(x.noise)
    S = 20 .* log10.(abs.(fft(w_noi)))
    lines!(ax_s, freqaxis(x.noise) ./ 1000, S .- (isnan(norm) ? maximum(S) : norm))
    xlims!(ax_s, xlims...)
    ylims!(ax_s, ylims...)

    w_sig = synthesize(x.signal)
    S = 20 .* log10.(abs.(fft(w_sig)))
    lines!(ax_s2, freqaxis(x.signal) ./ 1000, S .- (isnan(norm) ? maximum(S) : norm))
    xlims!(ax_s2, xlims...)
    ylims!(ax_s2, ylims...)
    
    ax_s.title = "Masker"
    ax_s2.title = "Signal"

    # Add labels
    ax_w.xlabel = "Time (s)"
    ax_w.ylabel = "Amplitude (Pa)"
    ax_s.xlabel = "Frequency (kHz)"
    ax_s.ylabel = "Level (dB re: max)"
    ax_s2.xlabel = "Frequency (kHz)"
    ax_s2.ylabel = "Level (dB re: max)"

    # Add title
    Label(fig[0, :], gettitle(x); fontsize=20.0)

    fig
end


# Commented out as we re-consider roved stimuli!
# """
#     RovedStimulus

# Wrapper for a regular stimulus that randomizes one or more stimulus parameters 

# Wrapper for typical subtypes of stimulus that implements the concept of "roving" or
# parameter randomization. RovedStimuli provide a `rand` method that returns a single sample
# from the (multivariate) distribution of the parameter(s) that is(are) randomized. It also
# provides an implementation of `synthesize` method that draws a parameter sample using `rand` 
# and then synthesizes a copy of the stimulus with those parameter values.
# """
# @with_kw struct RovedStimulus{S} <: RandomStimulus where {S <: AbstractStimulus}
#     rove_dist::Distribution
#     rove_params::Vector{Symbol}=Symbol[]
#     stimulus::S
# end

# Base.rand(s::RovedStimulus, args...; kwargs...) = rand(s.rove_dist, args...; kwargs...)

# RovedStimulus(s::AbstractStimulus, n::Int; kwargs...) = repeat([RovedStimulus(kwargs[:rove_dist], kwargs[:rove_params], s)], n)

# function Utilities.synthesize(s::RovedStimulus{S}) where {S <: AbstractStimulus}
#     # Get all field names and corresponding values for the underlying stimulus
#     keyvals = Dict([(k, getfield(s.stimulus, k)) for k in fieldnames(typeof(s.stimulus))])

#     # Draw a sample from the rove distribution and update relevant values
#     θ = rand(s)
#     @assert length(θ) == length(s.rove_params)
#     for (k, v) in zip(s.rove_params, θ)
#         keyvals[k] = v
#     end
#     # Synthesize a new stimulus copy
#     stimulus_copy = S(; keyvals...)

#     # Synthesize and return
#     synthesize(stimulus_copy)
# end

# function id(comp::RovedStimulus; accesses=nothing, connector="_", kwargs...)
#     id_main = savename(
#         string(typeof(comp)),
#         comp; 
#         accesses=accesses === nothing ? fieldnames(typeof(comp)) : accesses,
#         allowedtypes=(
#             Real, 
#             String, 
#             Symbol, 
#             Function,
#             Component,
#             Audiogram,
#         ), 
#         connector=connector,
#         kwargs...
#     )

#     names = join(string.(comp.rove_params), connector)
#     names = names * connector * "roved_over"
#     dist = string(comp.rove_dist)

#     return join([id_main, names, dist], connector)
# end
