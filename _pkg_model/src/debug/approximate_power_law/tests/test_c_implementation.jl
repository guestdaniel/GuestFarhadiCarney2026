using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

## Write function to side-by-side
function compare(stim=stim, stages=["syn", "hsr", "lsr"])
    # Simulate efferent-model response w/ true powerlaw
    @time resp_true = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=1)
    @time resp_appr = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=2)

    # Compute time axis 
    t = 0.0:(1/100e3):(length(resp_true["hsr"])/100e3 - 1/100e3)

    # Create LPF 
    responsetype = Lowpass(100.0; fs=100e3)
    designmethod = Butterworth(4)
    lpf = digitalfilter(responsetype, designmethod)

    # Set up plot
    fig = Figure()
    axs = [Axis(fig[i, j]) for i in 1:length(stages), j in 1:4]

    # Plot data
    for (idx, stage) in enumerate(stages)
        lines!(axs[idx, 1], t, resp_true[stage])
        lines!(axs[idx, 2], t, resp_appr[stage])
        lines!(axs[idx, 3], t, resp_appr[stage] .- resp_true[stage])
        lines!(axs[idx, 4], t, filt(lpf, resp_appr[stage]) .- filt(lpf, resp_true[stage]))
    end

    # Adjust limits
    xlims!.(axs, 0.0, 0.04)
    ylims!.(axs[1, 1:2], 0.0, 4500.0)
    ylims!.(axs[1, 3:4], -100.0, 100.0)
    ylims!.(axs[2, 1:2], 0.0, 1200.0)
    ylims!.(axs[2, 3:4], -40.0, 40.0)
    ylims!.(axs[3, 1:2], 0.0, 200.0)
    ylims!.(axs[3, 3:4], -10.0, 10.0)

    # Add labels
    Label(fig[0, 1], "True PLA"; tellwidth=false)
    Label(fig[0, 2], "Approximate PLA"; tellwidth=false)
    Label(fig[0, 3], "Difference\n(Approximate − true)"; tellwidth=false)
    Label(fig[0, 3], "LPF Difference\n(Approximate − true)"; tellwidth=false)

    # Neaten
    hidexdecorations!.(axs[1:2, :], grid=false, ticks=false, minorticks=false)
    hideydecorations!.(axs[:, 2], grid=false, ticks=false, minorticks=false)
    hideydecorations!.(axs[:, 4], grid=false, ticks=false, minorticks=false)

    fig, axs
end

## Use!
fig, axs = compare(scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 1.0, 100e3), 0.01, 100e3), 50.0))
fig
#xlims!.(axs, 0.01, 0.015)
#fig