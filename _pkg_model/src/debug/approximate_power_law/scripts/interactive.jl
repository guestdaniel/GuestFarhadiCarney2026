using Helios
using GLMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

freq = 8000.0
stim = vcat(scale_dbspl(cosine_ramp(pure_tone(freq, 0.0, 0.2, 100e3), 0.01, 100e3), 70.0), zeros(5000))
resp_true = sim_gfc2023_dict(stim, freq; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=1)
resp_appr = sim_gfc2023_dict(stim, freq; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=2)
expon = resp_true["expon"]
appro = adapt_ea_iir_parallel_sb(expon, 0.07669 * 2.5e-1, LogRange(5.7811e-4, 7.8475e3, 200); fs=100e3)[1]

fig = Figure()
ax = Axis(fig[1, 1])
t = 0.0:(1/100e3):(length(expon)/100e3 - 1/100e3)
lines!(ax, t, resp_true["sout1"]; color=:black)
lines!(ax, t, appro; color=:red)
#xlims!(ax, 0.0, 0.03)
fig