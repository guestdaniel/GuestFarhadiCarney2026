using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using Optim

# Simulate response to 1-kHz pure tone at 50 dB SPL
stim = vcat(scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.2, 100e3), 0.01, 100e3), 50.0), zeros(10000))
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=1)

# Plot
fig = Figure()
t = 0.0:(1/100e3):(length(stim)/100e3 - 1/100e3)
ax = Axis(fig[1, 1])
lines!(ax, t, resp["expon"])
ax = Axis(fig[1, 2])
lines!(ax, t, resp["sout1"])
ax = Axis(fig[1, 3])
lines!(ax, t, resp["sout2"])
fig