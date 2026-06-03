export PrecursorStimulus, AbstractPrecursorStimulus, 
       MixedStimulus, 
       PaddedStimulus,
       BinauralStimulus, AbstractBinauralStimulus,
       ipsi, contra


# Define AbstractPrecursorStimulus
#
# AbstractPrecursorStimulus concrete subtypes contain two stimuli: `stimulus` which is
# preceded in time by `precursor` with a specified interstimulus interval.
abstract type AbstractPrecursorStimulus <: AbstractCompoundStimulus end

# Overwrite `dur` to be defined as composite duration of precursor, ISI, and stimulus
dur(x::AbstractPrecursorStimulus) = dur(x.precursor) + x.dur_isi + dur(x.stimulus)

# Overwrite `idxswin` to point to target/post-precursor stimulus only
idxswin(x::AbstractPrecursorStimulus) = sampleat(onsets(x)[2], samprate(x)) .+ idxswin(x.stimulus)

# Overwrite `onsets` and ilk to point to precursor and stimulus onsets, etc.
onsets(x::AbstractPrecursorStimulus) = [0.0, dur(x.precursor) + x.dur_isi]
offsets(x::AbstractPrecursorStimulus) = [dur(x.precursor) - sampint(x), dur(x.precursor) + x.dur_isi + dur(x.stimulus) - sampint(x)]

# Define `synthesize` method 
synthesize(stimulus::AbstractPrecursorStimulus) = 
    vcat(
        synthesize(stimulus.precursor),
        zeros(samples(stimulus.dur_isi, stimulus.fs)),
        synthesize(stimulus.stimulus)
    )

# Define `parts`
parts(x::AbstractPrecursorStimulus) = Dict("precursor" => x.precursor, "probe" => x.stimulus)

# We provide a standard concrete subtype of AbstractPrecursorStimulus
"""
    PrecursorStimulus{A, B} <: AbstractCompoundStimulus

Compound type combining `precursor` and `stimulus` with silence of `dur_isi` 
"""
@with_kw struct PrecursorStimulus{A, B} <: AbstractPrecursorStimulus where {A <: AbstractStimulus, B <: AbstractStimulus}
    precursor::A
    stimulus::B
    fs::Float64=precursor.fs
    dur_isi::Float64=0.30
end

# Method for generating PrecursorStimulus from two stimuli
PrecursorStimulus(a, b; kwargs...) = PrecursorStimulus(; precursor=a, stimulus=b, kwargs...)

# Here we define another `viz` method but specifically for AbstractPrecursorStimulus
function viz(
    x::AbstractPrecursorStimulus; 
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
    w_pre = synthesize(x.precursor)
    S = 20 .* log10.(abs.(fft(w_pre)))
    lines!(ax_s, freqaxis(x.precursor) ./ 1000, S .- (isnan(norm) ? maximum(S) : norm))
    xlims!(ax_s, xlims...)
    ylims!(ax_s, ylims...)

    w_sig = synthesize(x.stimulus)
    S = 20 .* log10.(abs.(fft(w_sig)))
    lines!(ax_s2, freqaxis(x.stimulus) ./ 1000, S .- (isnan(norm) ? maximum(S) : norm))
    xlims!(ax_s2, xlims...)
    ylims!(ax_s2, ylims...)
    
    ax_s.title = "Precursor"
    ax_s2.title = "Probe"

    # Add labels
    ax_w.xlabel = "Time (s)"
    ax_w.ylabel = "Amplitude (Pa)"
    ax_s.xlabel = "Frequency (kHz)"
    ax_s.ylabel = "Level (dB re: max)"
    ax_s2.xlabel = "Frequency (kHz)"
    ax_s2.ylabel = "Level (dB re: max)"

    fig
end


struct PaddedStimulus{T} <: AbstractCompoundStimulus where {T <: AbstractStimulus}
    stim::T
    dur_pre::Float64
    dur_post::Float64
    fs::Float64
end

function PaddedStimulus(stim::AbstractStimulus, dur_pre::Float64, dur_post::Float64; fs=100e3)
    PaddedStimulus(stim, dur_pre, dur_post, fs)
end

function PaddedStimulus(stim::AbstractStimulus, dur_total::Float64; fs=100e3)
    dur_fringe = dur_total - dur(stim)
    @assert dur_fringe >= 0.0
    PaddedStimulus(stim, dur_fringe/2, dur_fringe/2, fs)
end

function synthesize(stim::PaddedStimulus) 
    pre = zeros(samples(stim.dur_pre, samprate(stim)))
    post = zeros(samples(stim.dur_post, samprate(stim)))
    return vcat(pre, synthesize(stim.stim), post)
end

samprate(stim::PaddedStimulus) = samprate(stim.stim)
dur(stim::PaddedStimulus) = stim.dur_pre + stim.dur_post + dur(stim.stim)
level(stim::PaddedStimulus) = level(stim.stim)
freq(stim::PaddedStimulus) = freq(stim.stim)

############################################################################################
# MixedStimulus
#
# MixedStimulus contains two stimuli: `a` and `b`, mixed evenly with a delay of `d` seconds
# imposed on `b`
#
# Here we define the type and implement several methods for it. 

"""
    MixedStimulus{A, B} <: AbstractStimulus

Compound type mixing stimuli `a` and `b` with delay `d` imposed on `b`
"""
@with_kw struct MixedStimulus{A, B} <: AbstractCompoundStimulus where {A <: AbstractStimulus, B <: AbstractStimulus}
    a::A
    b::B
    fs::Float64=a.fs
    d::Float64=0.0
end

# Method for generating MixedStimulus from two stimuli
MixedStimulus(a, b; d=0.0, kwargs...) = MixedStimulus(; a=a, b=b, d=d, kwargs...)

# Overwrite `dur` to indicate dur of `a`, which is expected to be longer 
dur(x::MixedStimulus) = dur(x.a)

# Overwrite `onsets` and related methods, pointing to onset of `b`
onsets(x::MixedStimulus) = [onset(x.a), onset(x.b) + x.d]
onset(x::MixedStimulus) = x.d
offset(x::MixedStimulus) = onset(x) + dur(x.b) - 1/samprate(x)

# `synthesize` method, returns mixed time-pressure waveform
synthesize(x::MixedStimulus) = mixat(synthesize(x.a), synthesize(x.b), x.d)

# Define `parts`
parts(x::MixedStimulus) = Dict("stim1" => x.a, "stim2" => x.b)


############################################################################################
# Define BinauralStimulus
abstract type AbstractBinauralStimulus <: AbstractCompoundStimulus end

"""
    BinauralStimulus{A, B} <: AbstractBinauralStimulus

Compound type containing two stimuli: `ipsi` in ipsilateral (first) channel, `contra` in contralateral (second) channel
"""
struct BinauralStimulus{A, B} <: AbstractBinauralStimulus where {A <: AbstractStimulus, B <: AbstractStimulus}
    ipsi::A
    contra::B
end

# Method for the same stimulus in both ears (becomes simple wrapper for diotic stimuli)
BinauralStimulus(x::AbstractStimulus) = BinauralStimulus(x, x)

# Overwrite `samperate` method to point to ipsi stimulus, and implement `dur`
samprate(x::AbstractBinauralStimulus) = samprate(ipsi(x))
dur(x::AbstractBinauralStimulus) = dur(ipsi(x))

# Methods `ipsi` and `contra` for accessing ipsi- and contra-side models
ipsi(x::AbstractBinauralStimulus) = x.ipsi
contra(x::AbstractBinauralStimulus) = x.contra

# `synthesize` method, returns tuple of left and right time-pressure waveforms
synthesize(stim::AbstractBinauralStimulus) = (synthesize(stim.ipsi), synthesize(stim.contra))

# `viz` method, plots ipsi on top and contra on bottom
function viz(stim::AbstractBinauralStimulus; fig=Figure(; size=(800, 500)))
    times = timeat.(extrema(idxswin(stim)), samprate(stim))
    viz(ipsi(stim); fig=fig, ax_w=Axis(fig[1, 1]), ax_s=Axis(fig[1, 2]), norm=50.0, show_title=false, times=times)
    viz(contra(stim); fig=fig, ax_w=Axis(fig[2, 1]), ax_s=Axis(fig[2, 2]), norm=50.0, show_title=false, times=times)
    Label(fig[1, 3], "Ipsilateral"; tellheight=false)
    Label(fig[2, 3], "Contralateral"; tellheight=false)
    fig
end

# Define `parts`
parts(x::AbstractBinauralStimulus) = Dict("ipsi" => ipsi(x), "contra" => contra(x))
