export PureTone, ComplexTone, SAMTone, SFMTone, ProfileAnalysisTone, SpectrotemporalRipple

"""
    PureTone <: AbstractStimulus

Pure tone at `freq` and phase `ϕ` scaled to overall `level` with `dur`
"""
@with_kw struct PureTone <: AbstractToneStimulus
    freq::Float64 = 1000.0
    ϕ::Float64 = 0π
    dur::Float64 = 1.0
    dur_ramp::Float64 = min(dur / 100, 0.01)
    level::Float64 = 50.0
    fs::Float64 = 100e3
end

PureTone(levels::AbstractVector; kwargs...) = [PureTone(; level=l, kwargs...) for l in levels]

function synthesize(x::PureTone)
    Parameters.@unpack freq, ϕ, dur, fs, dur_ramp, level = x
    scale_dbspl(cosine_ramp(pure_tone(freq, ϕ, dur, fs), dur_ramp, fs), level)
end

"""
    ComplexTone <: AbstractStimulus

Sum of pure tones, each individually variable in terms of level or other params
"""
@with_kw struct ComplexTone <: AbstractToneStimulus
    comps::Vector{PureTone} = PureTone[]
    level::Float64 = NaN
    dur::Float64 = 1.0
    dur_ramp::Float64 = min(dur / 100, 0.01)
    fs::Float64 = 100e3
end

function synthesize(x::ComplexTone)
    # Synthesize and combine all the individual components
    y = sum(map(synthesize, x.comps))

    # If overall level is NOT NaN, scale to overall level
    if !isnan(x.level)
        y = scale_dbspl(y, x.level)
    end
    return y
end

function ComplexTone(freqs, levels; level=NaN, fs=100e3, dur=1.0, dur_ramp=0.01, kwargs...)
    comps = [PureTone(freq=f, level=l; fs=fs, dur=dur, dur_ramp=dur_ramp, kwargs...) for (f, l) in zip(freqs, levels)]
    ComplexTone(; comps=comps, level=level, dur=dur, dur_ramp=dur_ramp, fs=fs)
end

freq(x::ComplexTone) = freq(first(x.comps))

"""
    SAMTone <: AbstractToneStimulus

SAM tone with carrier `freq` and mod `fm` and phase `ϕ` scaled to overall `level` with `dur`
"""
@with_kw struct SAMTone <: AbstractToneStimulus
    freq::Float64 = 1000.0
    fm::Float64 = 64.0
    ϕ::Float64 = 0π
    ϕ_m::Float64 = -π / 2
    m::Float64 = 1.0
    dur::Float64 = 1.0
    dur_ramp::Float64 = min(dur / 100, 0.01)
    level::Float64 = 50.0
    fs::Float64 = 100e3
end

function synthesize(x::SAMTone)
    @unpack_SAMTone x
    carrier = pure_tone(freq, ϕ, dur, fs)
    modulator = pure_tone(fm, ϕ_m, dur, fs)
    y = (1 .+ m .* modulator) .* carrier
    y = cosine_ramp(y, dur_ramp, fs)
    y = scale_dbspl(y, level)
    return y
end

modfreq(x::SAMTone) = x.fm

"""
    SFMTone <: AbstractToneStimulus

SFM tone with carrier `freq` and mod `fm` and phase `ϕ` scaled to overall `level` with `dur`. `depth` is specified in terms of Δf (Hz)
"""
@with_kw struct SFMTone <: AbstractToneStimulus
    freq::Float64 = 1000.0
    fm::Float64 = 64.0
    ϕ::Float64 = 0π
    ϕ_m::Float64 = -π / 2
    depth::Float64 = freq / 100
    dur::Float64 = 1.0
    dur_ramp::Float64 = min(dur / 100, 0.01)
    level::Float64 = 50.0
    fs::Float64 = 100e3
end

function synthesize(x::SFMTone)
    @unpack_SFMTone x

    # Manually synthesize using phase modulation
    t = timeaxis(x)
    Φ_base = 2π .* freq .* t .+ ϕ
    β = depth / fm
    Φ_mod = β .* sin.(2π .* fm .* t .+ ϕ_m)
    y = sin.(Φ_base .+ Φ_mod)

    # Ramp and scale for output
    y = cosine_ramp(y, dur_ramp, fs)
    y = scale_dbspl(y, level)
    return y
end

modfreq(x::SFMTone) = x.fm

"""
    ProfileAnalysisTone <: AbstractToneStimulus

Sum of pure tones in the style of Green profile analysis
"""
@with_kw struct ProfileAnalysisTone <: AbstractToneStimulus
    freqs::Vector{Float64}
    target_comp::Int = middle(eachindex(freqs))
    pedestal_level::Float64 = 50.0
    increment::Float64 = 0.0
    dur::Float64 = 0.3
    dur_ramp::Float64 = min(dur / 100, 0.01)
    fs::Float64 = 100e3
end

"""
    profile_analysis_tone(freqs::Vector, [target_comp::Int]; kwargs...)

Synthesizes a profile analysis tone composed of a deterministic set of components

Synthesized in the following way, per procedure used in Carney lab for so-called
`profile_analysis_iso` as of 09/22/2022:
- Background is synthesized (i.e., stimulus w/o target component), and the
  overall background level is set to the requested level. That is, the requested
  levels refers to the overall sound level of the background in dB SPL
- Target component is added in, with -infy dB SRS yielding a target component w/
  with the same amplitude as the background components and 0 dB yielding a target
  component with twice the amplitude of the background components

# Arguments
- `freqs::Tuple`: vector of frequenices to include in stimulus
- `target_comp`: which component should contain the increment (index into component_freqs, see code below)
- `fs=100e3`: sampling rate (Hz)
- `dur=0.10`: duration (s)
- `dur_ramp=0.01`: ramp duration (s)
- `pedestal_level=50.0`: overall sound level of background (dB SPL)
- `increment=0.0`: increment size in units of signal re: standard (dB)

# Returns
- `::Vector`: vector containing profile analysis tone
"""
function profile_analysis_tone(
    freqs,
    target_comp=Int(ceil(length(freqs) / 2));
    fs=100e3,
    dur=0.10,
    dur_ramp=0.01,
    pedestal_level=70.0,
    increment=0.0,
    phase_mode="fixed"
)
    # First, synthesize background, including component at target frequency
    background = map(freqs) do freq
        if phase_mode == "fixed"
            pure_tone(freq, 0.0, dur, fs)
        else
            pure_tone(freq, 2π * rand(), dur, fs)
        end
    end
    background = sum(background)

    # Given this stimulus, calculate required gain to achieve desired background
    # / pedestal level
    gain = pedestal_level - dbspl(background)

    # Synthesize target increment
    if phase_mode == "fixed"
        target = pure_tone(freqs[target_comp], 0.0, dur, fs) .* (10 .^ (increment ./ 20))
    else
        target = pure_tone(freqs[target_comp], 2π * rand(), dur, fs) .* (10 .^ (increment ./ 20))
    end

    # Add together, ramp, and scale
    stimulus = background .+ target
    stimulus = amplify(stimulus, gain)
    stimulus = cosine_ramp(stimulus, dur_ramp, fs)

    # Return
    return stimulus
end

function synthesize(x::ProfileAnalysisTone)
    @unpack_ProfileAnalysisTone x
    profile_analysis_tone(freqs, target_comp; fs=fs, dur=dur, dur_ramp=dur_ramp, pedestal_level=pedestal_level, increment=increment)
end



"""
    SpectrotemporalRipple <: AbstractToneStimulus

Spectrotemporal ripple stimulus.
"""
@with_kw struct SpectrotemporalRipple <: AbstractToneStimulus
    freq_low::Float64 = 250.0      # Hz
    freq_high::Float64 = 8000.0    # Hz
    cutoff_low::Float64 = 250.0    # Hz
    cutoff_high::Float64 = 8000.0  # Hz
    filter_order::Int = 8
    rate::Float64 = 1.0 # Hz
    density::Float64 = 1.0 # cyc/oct
    n_component::Int = 100
    level::Float64 = 65.0
    dur::Float64 = 1.0
    dur_ramp::Float64 = 0.01
    orientation::Int = 1 # 1 for upward, -1 for downward
    fs::Float64 = 100e3
    depth::Float64 = 0.0 # dB (0 == fully modulated)
end

function synthesize(x::SpectrotemporalRipple)
    @unpack_SpectrotemporalRipple x

    # Create filter
    responsetype = Bandpass(cutoff_low, cutoff_high; fs=fs)
    designmethod = Butterworth(div(filter_order, 4))
    f = digitalfilter(responsetype, designmethod)

    # Choose frequencies
    freqs = LogRange(freq_low, freq_high, n_component)

    # Time vector
    len_stim = round(Int, dur * fs)
    t = timevec(len_stim, fs)

    # Storage
    output_sig = zeros(len_stim)
    ripple_phase = 2π * rand()
    m = 10^(depth / 20)

    # Loop through components
    for (i, freq) in enumerate(freqs)
        # TODO: scrutinize this line --- still wack
        mod_phase = 2π * (orientation * rate .* t .+ density * log2(freq / freqs[1])) .+ 2π * ripple_phase
        modulator = 1 .+ m .* sin.(mod_phase)

        carrier_phase = 2π * rand()
        carrier = sin.(2π * freq .* t .+ carrier_phase)

        component = modulator .* carrier

        # Scale component to level
        component = scale_dbspl(component, level)

        output_sig .+= component
    end

    # Filter
    output_sig = filtfilt(f, output_sig)

    # Ramp
    output_sig = cosine_ramp(output_sig, dur_ramp, fs)

    return output_sig
end
