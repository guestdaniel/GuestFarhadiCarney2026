export ProbePureToneElicitor, ProbePureToneElicitor3

"""
    ProbePureToneElicitor <: AbstractCompoundStimulus
"""
struct ProbePureToneElicitor <: AbstractBinauralStimulus
    ipsi::PaddedStimulus
    contra::PaddedStimulus
    dur_probe::Float64
    dur_suppressor::Float64
    dur_wait::Float64
    dur_delay::Float64
    dur_post::Float64
    t_onset::Float64
end

function ProbePureToneElicitor(
    params_probe::Vector=[1e3, 20.0], 
    params_suppressor::Vector=[2e3, 60.0]; 
    dur_probe=0.1,
    dur_suppressor=0.1,
    dur_wait=0.0,
    dur_delay=-0.1,
    dur_post=0.0,
    fs=100e3,
)
    # Create probe
    freq, level = params_probe
    probe = PureTone(; freq=freq, level=level, dur=dur_probe, fs=fs)

    # Create suppressor, branching based on dur_delay
    t_onset = 0.0  # variable to record onset of probe 
    freq, level = params_suppressor
    if dur_delay < 0.0
        # In this case, the suppressor starts *before* the probe, so we need to make
        # the suppressor longer
        suppressor = PureTone(; freq=freq, level=level, dur=dur_suppressor+abs(dur_delay))

        # Then we need to append silence to the probe to match
        probe = PaddedStimulus(probe, abs(dur_delay), 0.0)

        # In this case, t_onset is
        t_onset = abs(dur_delay)
    else
        # In this case, the suppressor starts *after* the probe. 
        suppressor = PureTone(; freq=freq, level=level, dur=dur_suppressor)

        # We add delay to the suppressor
        suppressor = PaddedStimulus(suppressor, dur_delay, 0.0)

        # In this case, t_onset is just 0.0
        t_onset = 0.0
    end

    # If the probe and/or suppressor including silent gaps do not match in duration, we need 
    # to match them
    if dur(probe) > dur(suppressor)
        suppressor = PaddedStimulus(suppressor, 0.0, dur(probe)-dur(suppressor))
    elseif dur(probe) < dur(suppressor)
        probe = PaddedStimulus(probe, 0.0, dur(suppressor) - dur(probe))
    end

    # Finally, we need to add dur_post and dur_wait to both probe and suppressor
    probe = PaddedStimulus(probe, dur_wait, dur_post)
    suppressor = PaddedStimulus(suppressor, dur_wait, dur_post)
    t_onset = t_onset + dur_wait
    ProbePureToneElicitor(probe, suppressor, dur_probe, dur_suppressor, dur_wait, dur_delay, dur_post, t_onset)
end

# Overload idxswin for this stimulus to be post-onset window
function Utilities.idxswin(stim::ProbePureToneElicitor; dur_skip_onset=0.05)
    idx_onset = sampleat(stim.t_onset + dur_skip_onset, samprate(stim))
    idxs = idx_onset:(idx_onset + samples(stim.dur_probe - dur_skip_onset, samprate(stim)) - 1)
    @assert length(idxs)/samprate(stim) >= 0.05  # assert we have at least 50 ms of probe resp
    return idxs
end

# Overload other relevant bits
Utilities.level(stim::ProbePureToneElicitor) = level(ipsi(stim))
Utilities.freq(stim::ProbePureToneElicitor) = Utilities.freq(ipsi(stim))

"""
    ProbePureToneElicitor3 <: AbstractCompoundStimulus

Much like `ProbePureToneElicitor`, but with a different and generally more convenient 
parametrization. Here, we assume for simplicity that the probe always follows the elicitor
by `dur_delay` seconds. The sum of `dur_elicitor` must match the sum of `dur_probe` 
and `dur_delay`.
"""
struct ProbePureToneElicitor3 <: AbstractBinauralStimulus
    ipsi::PaddedStimulus
    contra::PureTone
    dur_probe::Float64
    dur_elicitor::Float64
    dur_delay::Float64
end

function ProbePureToneElicitor3(
    params_probe::Vector=[1e3, 20.0], 
    params_elicitor::Vector=[2e3, 60.0]; 
    dur_probe=0.1,
    dur_elicitor=1.0,
    dur_delay=dur_elicitor - dur_probe,
    dur_ramp=0.005, 
    fs=100e3,
)
    # Create probe
    freq, level = params_probe
    probe = PureTone(; freq=freq, level=level, dur=dur_probe, fs=fs, dur_ramp=dur_ramp)
    probe = PaddedStimulus(probe, dur_delay, 0.0, fs)

    # Create elicitor
    freq, level = params_elicitor
    elicitor = PureTone(; freq=freq, level=level, dur=dur_elicitor, fs=fs, dur_ramp=dur_ramp)

    # Finally, we need to add dur_post and dur_wait to both probe and suppressor
    ProbePureToneElicitor3(probe, elicitor, dur_probe, dur_elicitor, dur_delay)
end

# Overload idxswin for this stimulus to be post-onset window
function Utilities.idxswin(stim::ProbePureToneElicitor3; dur_skip_onset=0.05)
    idx_onset = sampleat(stim.dur_delay + dur_skip_onset, samprate(stim))
    idxs = idx_onset:(idx_onset + samples(stim.dur_probe - dur_skip_onset, samprate(stim)) - 1)
    @assert length(idxs)/samprate(stim) >= 0.05  # assert we have at least 50 ms of probe resp
    return idxs
end

# Overload other relevant bits
Utilities.level(stim::ProbePureToneElicitor3) = level(ipsi(stim))
Utilities.freq(stim::ProbePureToneElicitor3) = Utilities.freq(ipsi(stim))