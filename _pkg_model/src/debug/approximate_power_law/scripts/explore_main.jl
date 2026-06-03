using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

## Handle prep work (synthesize inputs, set parameter values, etc.)
# Set params
fs = 100e3

# Synthesize 1-kHz pure tone at 50 dB SPL
stim = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 1.0, 100e3), 0.01, 100e3), 50.0)

# Simulate efferent-model response
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0)

# Extract and downsample output waveform from exponential adaptation stage (input to PLA)
resp_old = sim_orig_dict(stim, 1000.0)
x = resample(resp["expon"], fs/100e3)
x = vcat(x, zeros(5000))
t = 0.0:(1/fs):(length(x)/fs - 1/fs)

# Set param values
α = 2.5e-1
τ_short=5e-3
τ_long=1e4

# Simulate true power-law adaptation
@time y_pla, _ = adapt_pla(x, α, 5e-3; fs=fs)

# Define f and g
function f(θ)
    y_ea, _, _ = adapt_ea_iir_parallel(x, θ[1] * α, LogRange(exp(θ[2]), exp(θ[3]), 125); fs=fs)
    rms(y_pla .- y_ea)
end

x_init = [1.0, log(τ_short), log(τ_long)]
results = optimize(f, x_init)
x̂ = Optim.minimizer(results)

ŷ, _, _ = adapt_ea_iir_parallel(x, x̂[1] * α, LogRange(exp(x̂[2]), exp(x̂[3]), 125); fs=fs)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, t, y_pla; color=:black)
lines!(ax, t, ŷ; color=:pink)
xlims!(0.0, 0.1)
fig