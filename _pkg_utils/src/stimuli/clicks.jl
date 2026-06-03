export Click, ClickSequence

"""
    Click <: AbstractStimulus

Rectangular click stimulus with specified duration and peak-equivalent SPL.

# Fields
- `level::Float64`: Peak-equivalent SPL (dB SPL). 
- `click_dur::Float64`: Duration of the click pulse (s). Default 100e-6.
- `dur::Float64`: Total duration of the stimulus buffer (s). Default 0.1.
- `onset::Float64`: Onset time of the click (s). Default 0.0.
- `fs::Float64`: Sampling rate (Hz). Default 100e3.
- `polarity::Float64`: Polarity of the click (1.0 or -1.0). Default 1.0.
"""
@with_kw struct Click <: AbstractStimulus
    level::Float64 = 100.0
    click_dur::Float64 = 100e-6
    dur::Float64 = 0.1
    onset::Float64 = 0.0
    fs::Float64 = 100e3
    polarity::Float64 = 1.0
end

function synthesize(x::Click)
    Parameters.@unpack level, click_dur, dur, onset, fs, polarity = x

    # Calculate peak amplitude from peSPL
    # peSPL = 20 * log10(peak_amp / (sqrt(2) * 20e-6))
    # peak_amp = sqrt(2) * 20e-6 * 10^(peSPL / 20)
    ref_pressure = 20e-6
    amp = sqrt(2) * ref_pressure * 10^(level / 20)

    # Create silence buffer
    w = zeros(samples(dur, fs))

    # Add click
    i_start = sampleat(onset, fs)
    i_end = sampleat(onset + click_dur, fs)

    # Ensure indices are valid
    i_start = max(1, i_start)
    i_end = min(length(w), i_end)

    if i_start <= i_end
        w[i_start:i_end] .= polarity * amp
    end

    return w
end

onset(x::Click) = x.onset
offset(x::Click) = x.onset + x.click_dur

"""
    ClickSequence <: AbstractStimulus

Sequence of clicks with increasing levels.

# Fields
- `n_clicks::Int`: Number of clicks in sequence (1-5). Default 5.
- `click_dur::Float64`: Duration of each click (s). Default 100e-6.
- `ici::Float64`: Inter-click interval (s). Default 1e-3.
- `start_level::Float64`: Level of first click (dB peSPL). Default 45.0.
- `step_level::Float64`: Level increment per click (dB). Default 15.0.
- `dur::Float64`: Total duration (s). Default 0.1.
- `onset::Float64`: Onset time of first click (s). Default 0.0.
- `fs::Float64`: Sampling rate (Hz). Default 100e3.
- `polarity::Float64`: Polarity of clicks. Default 1.0.
"""
@with_kw struct ClickSequence <: AbstractStimulus
    n_clicks::Int = 5
    click_dur::Float64 = 100e-6
    ici::Float64 = 1e-3
    start_level::Float64 = 45.0
    step_level::Float64 = 15.0
    dur::Float64 = 0.1
    onset::Float64 = 0.0
    fs::Float64 = 100e3
    polarity::Float64 = 1.0
end

function synthesize(x::ClickSequence)
    Parameters.@unpack n_clicks, click_dur, ici, start_level, step_level, dur, onset, fs, polarity = x

    # Create silence buffer
    w = zeros(samples(dur, fs))

    for i in 1:n_clicks
        # Calculate level and amplitude
        level = start_level + (i - 1) * step_level
        ref_pressure = 20e-6
        amp = sqrt(2) * ref_pressure * 10^(level / 20)

        # Calculate timing
        t_onset = onset + (i - 1) * ici
        i_start = sampleat(t_onset, fs)
        i_end = sampleat(t_onset + click_dur, fs)

        # Ensure indices are valid
        i_start = max(1, i_start)
        i_end = min(length(w), i_end)

        if i_start <= i_end
            w[i_start:i_end] .= polarity * amp
        end
    end

    return w
end

onsets(x::ClickSequence) = x.onset .+ (0:(x.n_clicks-1)) .* x.ici
offsets(x::ClickSequence) = onsets(x) .+ x.click_dur
onset(x::ClickSequence) = x.onset
offset(x::ClickSequence) = onsets(x)[end] + x.click_dur
