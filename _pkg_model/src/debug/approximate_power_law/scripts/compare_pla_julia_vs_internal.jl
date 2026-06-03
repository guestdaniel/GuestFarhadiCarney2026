using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
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
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, fs=100e3)

# Write functions for plotting
function plotcurve!(ax, x, y, title=""; kwargs...)
    xlims = (0.0, 0.01)
    ax.title = title
    lines!(ax, x, y; kwargs...)
    xlims!(xlims...)
end

# Write function for analysis
function analyze(resp)
    # Extract and output waveforms 
    expon = resp["expon"]
    y_external = adapt_pla_clike(expon, alpha1, beta1; fs=100e3)
    y_ccall = adapt_pla_c(expon, alpha1, beta1; fs=100e3)
    y_internal = resp["sout1"]
    t = 0.0:(1/100e3):(length(y_external)/100e3 - 1/100e3)

    # Start plot
    fig = Figure(; resolution=(800, 1500))

    # Plot Julia vs internal PLA
    plotcurve!(Axis(fig[1, 1]), t, stim, "Stimulus"; color=:black)
    plotcurve!(Axis(fig[2, 1]), t, expon, "Exponential-adaptation output (expon)"; color=:black)
    plotcurve!(Axis(fig[3, 1]), t, y_internal, "Powerlaw-adaptation output (sout1)"; color=:black)
    plotcurve!(Axis(fig[4, 1]), t, y_external, "Julia PLA (native)"; color=:black)
    plotcurve!(Axis(fig[5, 1]), t, y_ccall, "Julia PLA (ccall)"; color=:black)
    fig

    ax = Axis(fig[4, 2])
    plotcurve!(ax, t, y_internal; color=:orange)
    plotcurve!(ax, t, y_external, "Julia PLA (native)"; color=:black)
    ax = Axis(fig[5, 2])
    plotcurve!(ax, t, y_internal; color=:orange)
    plotcurve!(ax, t, y_ccall, "Julia PLA (ccall)"; color=:black)
    fig
end

analyze(resp)