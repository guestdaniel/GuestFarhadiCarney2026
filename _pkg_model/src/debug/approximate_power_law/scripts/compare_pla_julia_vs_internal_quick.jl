using Helios
using AuditorySignalUtils
using DSP
using Printf

## Compare PLA applied to exponential-adaptation output (expon) inside/outside of C code
# Set params 
alpha1 = 2.5e-6*100e3; 
beta1  = 5e-4; 

# Synthesize 1-kHz pure tone at 50 dB SPL
stim = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.2, 100e3), 0.01, 100e3), 50.0)

# Zero-pad
stim = vcat(stim, zeros(Int(round(0.1*100e3))))

# Compute model response
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, fs=100e3, dur_pad_left=0.0, dur_pad_right=0.0, clip_left=false, clip_right=false);
expon = resp["expon"];
y_ccall = adapt_pla_c(expon, alpha1, beta1; fs=100e3);

