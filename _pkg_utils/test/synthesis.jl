# Check to make sure that threshold-equalizing noise is the correct level
# If we synthesize threshold equalizing noise where the passband is exactly equal to that
# suggested by the value of the ERB at 1 kHz, then the level of the output waveform should
# nearly exactly match the level of the noise (which is specified in terms of the level in
# the 1 kHz ERB).
@test begin
    x = TENoise(; freq_low=0.935e3, freq_high=1.0681e3, dur_ramp=0.01, dur=5.0) 
    isapprox(level(x), dbspl(synthesize(x)); atol=0.1)
end

# Check to make sure that threshold-equalizing noise level increases with increasing passband width
# If we widen threshold-equalizing noise bandwidth, the level should increase, roughly as
# expected for a spectrally flat noise over a limited frequency range
@test begin
    x1 = TENoise(; freq_low=1e3, freq_high=2e3, dur_ramp=0.001, dur=5.0) 
    x2 = TENoise(; freq_low=1e3, freq_high=4e3, dur_ramp=0.001, dur=5.0) 
    isapprox(dbspl(synthesize(x2)), dbspl(synthesize(x1)) + 3; atol=1.0)
end