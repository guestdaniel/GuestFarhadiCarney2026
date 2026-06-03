using Distributed
addprocs(4)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Helios

using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

# Simulate efferent-model response w/ true powerlaw
@everywhere function test(dur, powerlaw_mode)
    @elapsed resp_true = sim_gfc2023_dict(zeros(Int(round(dur*100e3))), 1000.0; powerlaw_mode=powerlaw_mode)
end

# Collect data
durs = LogRange(0.01, 2.0, 20)
dur_true = pmap(x -> test(x, 1), durs)
dur_appr = pmap(x -> test(x, 2), durs)

# Create plot
fig = Figure(; resolution=(1100, 400))
ax = Axis(
    fig[1, 1], 
    xscale=log10, 
    xminorticksvisible=true,
    xminorticks=IntervalsBetween(9),
    xminorgridvisible=true,
    xticks = [0.01, 0.1, 1.0, 10.0],
    yscale=log10,
    yminorticksvisible=true,
    yminorticks=IntervalsBetween(9),
    yminorgridvisible=true,
    yticks = [0.01, 0.1, 1.0, 10.0],
)
lines!(ax, durs, dur_true)
lines!(ax, durs, dur_appr)
xlims!(ax, 0.009, 11.0)
ax.xlabel = "Stimulus time (s)"
ax.ylabel = "Compute time (s)"

ax = Axis(
    fig[1, 2], 
    xscale=log10, 
    xminorticksvisible=true,
    xminorticks=IntervalsBetween(9),
    xminorgridvisible=true,
    xticks = [0.01, 0.1, 1.0, 10.0],
    yscale=log10,
    yminorticksvisible=true,
    yminorticks=IntervalsBetween(9),
    yminorgridvisible=true,
    yticks = [0.01, 0.1, 1.0, 10.0],
)
lines!(ax, durs, dur_true ./ durs)
lines!(ax, durs, dur_appr ./ durs)
xlims!(ax, 0.009, 11.0)
ylims!(0.01, 12.0)
ax.xlabel = "Stimulus time (s)"
ax.ylabel = "Compute/stimulus time ratio"

ax = Axis(
    fig[1, 3], 
    xscale=log10, 
    xminorticksvisible=true,
    xminorticks=IntervalsBetween(9),
    xminorgridvisible=true,
    xticks = [0.01, 0.1, 1.0, 10.0],
    yscale=log10,
    yminorticksvisible=true,
    yminorticks=IntervalsBetween(9),
    yminorgridvisible=true,
    yticks = [0.01, 0.1, 1.0, 10.0],
)
lines!(ax, durs, dur_true ./ dur_appr; color=:black)
xlims!(ax, 0.009, 11.0)
ylims!(0.1, 100.0)
ax.xlabel = "Stimulus time (s)"
ax.ylabel = "Performance boost (ratio)"

fig