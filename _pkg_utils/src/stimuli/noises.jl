export GaussianNoise, GaussianNoiseSL, BandstopNoise, BandstopNoiseSL, SAMNoise, SAMNoiseSL, TENoise, NotchedNoiseSL, LFN, SAMLFN

"""
    GaussianNoise <: AbstractNoiseStimulus

Bandlimited Gaussian noise filtered in the spectral domain
"""
@with_kw struct GaussianNoise <: AbstractNoiseStimulus
    freq_low::Float64=20.0
    freq_high::Float64=10e3
    dur::Float64=0.5
    dur_ramp::Float64=0.01
    level::Float64=30.0 
    fs::Float64=100e3
end

function synthesize(x::GaussianNoise)
    # Unpack parameters
    Parameters.@unpack freq_low, freq_high, dur, fs, dur_ramp, level = x

    # Synthesize broadband Gaussian noise carrier
    w = randn(samples(dur, fs))

    # Filter in spectral domain
    w = specfilt_bp(w, freq_low, freq_high, fs)

    # Scale and ramp
    scale_dbspl!(w, level)
    cosine_ramp!(w, dur_ramp, fs)
    return w
end

"""
    GaussianNoiseSL <: AbstractNoiseStimulus

Bandlimited Gaussian noise filtered in the spectral domain and specified in spectrum level
"""
@with_kw struct GaussianNoiseSL <: AbstractNoiseStimulus
    freq_low::Float64=20.0
    freq_high::Float64=10e3
    dur::Float64=0.5
    dur_ramp::Float64=0.01
    level::Float64=30.0 
    fs::Float64=100e3
end

function synthesize(x::GaussianNoiseSL)
    # Unpack parameters
    Parameters.@unpack freq_low, freq_high, dur, fs, dur_ramp, level = x

    # Synthesize broadband Gaussian noise carrier
    w = randn(samples(dur, fs))

    # Filter in spectral domain
    w = specfilt_bp(w, freq_low, freq_high, fs)

    # Scale and ramp
    scale_dbspl!(w, sl_to_ol(level, freq_high-freq_low))
    cosine_ramp!(w, dur_ramp, fs)
end

"""
    LFN <: AbstractNoiseStimulus

Bandlimited low-fluctuation Gaussian noise filtered in the spectral domain
"""
@with_kw struct LFN <: AbstractNoiseStimulus
    freq_low::Float64=5000.0
    freq_high::Float64=5100.0
    dur::Float64=0.5
    dur_ramp::Float64=0.01
    level::Float64=30.0 
    fs::Float64=100e3
    N_iter::Int=10
end

function synthesize(x::LFN)
    # Unpack parameters
    Parameters.@unpack freq_low, freq_high, dur, fs, dur_ramp, level = x

    # Synthesize broadband Gaussian noise carrier
    w = randn(samples(dur, fs))

    # Filter in spectral domain
    w = specfilt_bp(w, freq_low, freq_high, fs)

    # Iteratively flatten the fluctuations by dividing by inst. Hilbert envelop and repeat filtering
    for _ in 1:x.N_iter
        env = abs.(hilbert(w))
        w .= w ./ env
        w = specfilt_bp(w, freq_low, freq_high, fs)
    end

    # Scale to overall level and ramp
    scale_dbspl!(w, level)
    cosine_ramp!(w, dur_ramp, fs)

    return w
end

"""
    BandstopNoise <: AbstractNoiseStimulus
"""
@with_kw struct BandstopNoise <: AbstractNoiseStimulus
    dur::Float64=0.5
    dur_ramp::Float64=0.01
    level::Float64=30.0 
    fs::Float64=100e3
    freq_low::Float64=20.0     # lower edge of the overall noise band
    freq_high::Float64=20e3    # upper edge of the overall noise band
    notch_low::Float64=900.0   # lower edge of the notch
    notch_high::Float64=1100.0 # upper edge of the notch
end

function synthesize(x::BandstopNoise)
    # Unpack parameters
    Parameters.@unpack freq_low, freq_high, dur, fs, dur_ramp, level, notch_low, notch_high = x

    # Synthesize broadband Gaussian noise carrier
    w = randn(samples(dur, fs))

    # Filter in spectral domain to create overall bandlimits
    w = specfilt_bp(w, freq_low, freq_high, fs)

    # Filter in spectral domain again to remove notch
    w = specfilt_br(w, notch_low, notch_high, fs)

    # Scale and ramp
    scale_dbspl!(w, level)
    cosine_ramp!(w, dur_ramp, fs)
    return w
end

"""
    BandstopNoiseSL <: AbstractNoiseStimulus
"""
@with_kw struct BandstopNoiseSL <: AbstractNoiseStimulus
    dur::Float64=0.5
    dur_ramp::Float64=0.01
    level::Float64=30.0 
    fs::Float64=100e3
    freq_low::Float64=20.0     # lower edge of the overall noise band
    freq_high::Float64=20e3    # upper edge of the overall noise band
    notch_low::Float64=900.0   # lower edge of the notch
    notch_high::Float64=1100.0 # upper edge of the notch
end

function synthesize(x::BandstopNoiseSL)
    # Unpack parameters
    Parameters.@unpack freq_low, freq_high, dur, fs, dur_ramp, level, notch_low, notch_high = x

    # Synthesize broadband Gaussian noise carrier
    w = randn(samples(dur, fs))

    # Filter in spectral domain
    w = specfilt_bp(w, freq_low, freq_high, fs)

    # Scale to correct level before doing notch reject; the logic is that spectrum level
    # can be most easily calculated from the overall bandwidth before the notch
    scale_dbspl!(w, sl_to_ol(level, freq_high-freq_low))

    # Filter in spectral domain again to remove notch
    w = specfilt_br(w, notch_low, notch_high, fs)

    # Ramp final signal
    cosine_ramp!(w, dur_ramp, fs)
end

"""
    TENoise <: AbstractNoiseStimulus

Threshold-equalizing noise, adapated from original code from lab of A. J. Oxenham.

Threshold-equalizing noise is originally described in:

Moore, B. C. J., Huss, M., Vickers, D. A., Glasberg, B. R., and Alcántra, J. I. (2000). “A
test for the diagnosis of dead regions in the cochlea,” British Journal of Audiology, 34,
205–224. doi:10.3109/03005364000000131

Threshold-equalizing noise (TEN) is a noise spectrally shaped so as to acheive uniform
masked thresholds for a pure tone as a function of frequency from 125–15000 Hz. TEN is 
produced under the assumption that the power of a signal at masked threshold is:

    P₀ = N₀ × K × ERB

where N₀ is the noise power spectral density, K is the signal-to-noise ratio at the output 
of the auditory filter required for masked threshold at the given frequency (this varies
with frequency), and ERB is the equivalent rectangular bandwidth of the auditory filter
at the given frequency.

Note that TEN specifies `level` in terms of level in the ERB centered at 1 kHz, and NOT in
terms of overall level. Thus, unlike for most other stimuli, `level(x)` does not return the
same number as `dbspl(synthesize(x))`. The underlying synthesis routine is provided by
`te_noise` in the AuditorySignalUtils package.
"""
@with_kw struct TENoise <: AbstractNoiseStimulus
    dur::Float64=0.1
    dur_ramp::Float64=0.01
    fs::Float64=100e3
    freq_low::Float64=125.0
    freq_high::Float64=15e3
    level::Float64=20.0  # note: level in 1 kHz ERB, NOT overall level
end

function synthesize(x::TENoise)
    # Unpack fields of TENoise object automagically
    @unpack_TENoise x

    te_noise(; fs=fs, dur=dur, level=level, dur_ramp=dur_ramp, freq_low=freq_low, freq_high=freq_high)
end

"""
    SAMNoise <: AbstractNoiseStimulus

Sinusoidally amplitude-modulated noise. 

Noise can be any type of bandlimited noise, kwargs are captured and sent to the appropriate
noise type based on the constructor defined below. In all cases, overall level is set at the
end following modulation. Note that a downside of this class currently is that there is 
limited control over how ramps interact with AM envelope.
"""
@with_kw struct SAMNoise{N} <: AbstractNoiseStimulus
    fm::Float64=10.0
    ϕm::Float64=π/2
    m::Float64=1.0
    carrier::N
end

function SAMNoise(carriertype::Type{T}; fm=10.0, ϕm=-π/2, m=1.0, kwargs...) where T<:AbstractNoiseStimulus
    carrier = T(; kwargs...)
    return SAMNoise{T}(;
        fm=fm,
        ϕm=ϕm,
        m=m,
        carrier=carrier,
    )
end

function synthesize(x::SAMNoise)
    # Synthesize carrier
    w = synthesize(x.carrier)

    # Synthesize modulator
    if x.fm == 0.0
        modulator = zeros(length(w))
    else
        modulator = pure_tone(x.fm, x.ϕm, x.carrier.dur, x.carrier.fs)
    end

    # Combine carrier and modulator and scale
    w = (1.0 .+ x.m .* modulator) .* w
    scale_dbspl!(w, x.carrier.level)
    return w
end

samprate(x::SAMNoise) = samprate(x.carrier)
dur(x::SAMNoise) = dur(x.carrier)