export pt, sam, moc_nonlinearity

# Utility function to synthesize pure tone quickly
pt(f=1000.0, l=50.0, dur=0.2, fs=100e3) = scale_dbspl(cosine_ramp(pure_tone(f, 0.0, dur, fs), 0.01, fs), l)

# Utility function to synthesize SAM tone quickly
function sam(f=1000.0, fm=1.0, d=0.0, l=50.0, dur=0.2, fs=100e3)
    m = 10^(d/10)
    carrier = pure_tone(f, 0.0, dur, fs)
    modulator = pure_tone(fm, 0.0, dur, fs)
    sam = (1.0 .+ m .* modulator) .* carrier
    sam = scale_dbspl(cosine_ramp(sam, 0.01, fs), l)
end

# Rational nonlinearity used in the MOC stage
function moc_nonlinearity(x, β=0.01, offset=0.0, maxrate=100.0, minrate=0.0) 
    if x < offset
        return maxrate
    else
        return ((maxrate-minrate) * 1.0/(1.0 + (β*(x-offset))^2)) + minrate
    end
end