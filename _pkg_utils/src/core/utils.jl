export spectrum, synchronize_cache, middle

"""
    spectrum(x, fs)
    spectrum(x, fs[; pad_factor=1])

Calculates spectrum of `x` given sampling rate of `fs`

Calculates spectrum of `x` using FFTW given a sampling rate of `fs`. Returns a tuple
containing a frequency-axis vector and the (two-sided) log-power spectrum. Has optional keyword argument
`pad_factor`. `pad_factor` of 1 yields just `fft(x)`, but for larger pad factors we zero-pad
by that factor before calculating FFT to implement ideal interpolation.
"""
function spectrum(x::Vector{T}, fs::T; pad_factor::Int=1) where {T <: Real}
    # Calculate frequency axis based on requested pad factor
    f = collect(LinRange(0.0, fs, length(x)*pad_factor))

    # Calculate window scaling factor
    S = length(x)

    # Calculate FFT of zero-padded input
    X = fft(vcat(x, zeros(length(x)*(pad_factor-1))))

    # Calculate power spectral density
    # Note that we take the magnitude spectrum, normalize it by 20 μPa, square it, and then
    # correct by a factor of 2/(fs*S) to arrive at the output power spectral density in 
    # units of dB SPL / Hz
#    psd = (1/fs) .* 2 .* (1/S) .* abs.(X ./ 20e-6).^2
    psd = 2 .* (1/S^2) .* abs.(X ./ 20e-6).^2
    # TODO Figure out --- does 2 go inside parentheses or outside? (i think it belongs inside!)

    return f, 10 .* log10.(psd)
end

"""
    synchronize_cache()
"""
function synchronize_cache()
    # Pull files from server first
    run(`wsl rsync -rv dguest2@bluehive.circ.rochester.edu:/scratch/dguest2/cl_cache /mnt/c/home/daniel`)
    # Send same files to server next
#    run(`wsl rsync -rv /mnt/c/home/daniel/cl_cache dguest2@bluehive.circ.rochester.edu:/scratch/dguest2`)
end

middle(x) = x[Int(ceil(length(x)/2))]
idxmiddle(x) = Int(ceil(length(x)/2))