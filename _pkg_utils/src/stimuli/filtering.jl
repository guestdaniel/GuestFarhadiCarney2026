export specfilt_bp, specfilt_br

"""
    specfilt_bp(x, lo, hi, fs=100e3)

Band-pass filter in the spectral domain between lo and hi
"""
function specfilt_bp(x, lo, hi, fs=100e3)
    # Determine fft frequencies
    f = fftfreq(length(x), fs)
    fabs = abs.(f)

    # Compute FFT
    X = fft(x)

    # Determine the frequency bins that fall within our band-reject area
    reject_bins = findall(.!(fabs .>= lo .&& fabs .<= hi))
   
    # Zero out those frequency components
    X[reject_bins] .= 0.0 + 0.0im

    # Resynthesize the time-domain signal
    return real.(ifft(X))
end

"""
    specfilt_br(x, lo, hi, fs=100e3)

Band-reject filter in the spectral domain between lo and hi
"""
function specfilt_br(x, lo, hi, fs=100e3)
    # Determine fft frequencies
    f = fftfreq(length(x), fs)
    fabs = abs.(f)

    # Compute FFT
    X = fft(x)

    # Determine the frequency bins that fall within our band-reject area
    band_reject_bins = findall(fabs .>= lo .&& fabs .<= hi)
   
    # Zero out those frequency components
    X[band_reject_bins] .= 0.0 + 0.0im

    # Resynthesize the time-domain signal
    return real.(ifft(X))
end