export gammatone_coef, filter_gammatone, filterbank_gammatone

"""
    gammatone_coef(f_c::Real, fs::Real, order::Integer=4)

Calculates the coefficients for a digital gammatone filter à la Patterson et al. (1987).

The gammatone function is:
    g(t) = t^(n-1) * exp(-2πb*t) * cos(2πf_c*t + ϕ)

This code produces coefficients for a cascade of four 2nd order IIR filters that approximate
convolution with the gammatone function. The original code for this implementation is in
MATLAB and available at: https://www.mathworks.com/matlabcentral/fileexchange/23053-gammatone-based-auditory-spectrograms

Patterson, R. D., Nimmo-Smith, I., Holdsworth, J., & Rice, P. (1987, December). An efficient
auditory filterbank based on the gammatone function. In a meeting of the IOC Speech Group on
Auditory Modelling at RSRE (Vol. 2, No. 7).

Slaney, M. (1993). An efficient implementation of the Patterson-Holdsworth auditory filter
bank. Apple Computer, Perception Group, Tech. Rep, 35(8).

Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory filter shapes from
notched-noise data. Hearing Research, 47(1-2), 103-138.
"""
function gammatone_coef(cf, fs; bw_factor=1.0)
    # Calculate quantites we need later
    T = 1/fs                                # Sampling interval
    bw_erb = 24.7 * (4.37 * cf / 1000 + 1)  # Equivalent rectangular bandwidth in Hz, Glasberg and Moore (1990)
    bw = 2π * bw_factor * 1.019 * bw_erb    # Bandwidth scaled as in Patterson et al. (1987) and specified param and multiplied by 2π

    # Create the coefficients for the gammatone filters
    A0 = T;    # A0 and A2 are constant 
    A2 = 0.0;
    A11 = -(T * cos(2*cf*pi*T) + T * sqrt(3 + 2^(1.5)) * sin(2*cf*pi*T)) ./ exp(bw*T);
    A12 = -(T * cos(2*cf*pi*T) - T * sqrt(3 + 2^(1.5)) * sin(2*cf*pi*T)) ./ exp(bw*T);
    A13 = -(T * cos(2*cf*pi*T) + T * sqrt(3 - 2^(1.5)) * sin(2*cf*pi*T)) ./ exp(bw*T);
    A14 = -(T * cos(2*cf*pi*T) - T * sqrt(3 - 2^(1.5)) * sin(2*cf*pi*T)) ./ exp(bw*T);

    B0 = 1;  # B0, B2, and B3 are constants
    B1 = -2 * cos(2*cf*pi*T) ./ exp(bw*T);
    B2 = exp(-2*bw*T);

    # Create gains for the gammatone filters
    gain = abs((-2*exp(4im*cf*pi*T)*T + 2*exp(-(bw*T) + 2*1im*cf*pi*T) * T * (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* sin(2*cf*pi*T))) *
           (-2*exp(4im*cf*pi*T)*T + 2*exp(-(bw*T) + 2*1im*cf*pi*T) * T * (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * sin(2*cf*pi*T))) *
           (-2*exp(4im*cf*pi*T)*T + 2*exp(-(bw*T) + 2*1im*cf*pi*T) * T * (cos(2*cf*pi*T) - sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) *
           (-2*exp(4im*cf*pi*T)*T + 2*exp(-(bw*T) + 2*1im*cf*pi*T) * T * (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) /
           (-2 / exp(2*bw*T) - 2*exp(4im*cf*pi*T) +  2*(1 + exp(4im*cf*pi*T)) / exp(bw*T))^4)

    # Return filter coefficients
    (a1, b1) = ([A0, A11, A2] ./ gain, [B0, B1, B2])
    (a2, b2) = ([A0, A12, A2], [B0, B1, B2])
    (a3, b3) = ([A0, A13, A2], [B0, B1, B2])
    (a4, b4) = ([A0, A14, A2], [B0, B1, B2])
    return (a1, b1), (a2, b2), (a3, b3), (a4, b4)
end

function filter_gammatone(y, cf, fs; kwargs...)
    filters = gammatone_coef(cf, fs; kwargs...)
    for f in filters
        y = filt(f..., y)
    end
    return y
end

filterbank_gammatone(y, cfs, fs; kwargs...) = [filter_gammatone(y, cf, fs; kwargs...) for cf in cfs]