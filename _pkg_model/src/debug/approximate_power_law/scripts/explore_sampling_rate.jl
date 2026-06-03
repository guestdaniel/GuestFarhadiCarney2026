using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

## Handle prep work (synthesize inputs, set parameter values, etc.)
# Synthesize 1-kHz pure tone at 50 dB SPL
stim = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 1.0, 100e3), 0.01, 100e3), 50.0)

# Simulate efferent-model response
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0)

# Extract and downsample output waveform from exponential adaptation stage (input to PLA)
x = resp["expon"]
x = vcat(x, zeros(5000))
t = 0.0:(1/100e3):(length(x)/100e3 - 1/100e3)

## Figure #2: Compare adaptation integrals and outputs for IIR parallel vs real PLA
function fig2(x=x, scale=1.0, τ_short=5e-3, τ_long=1e4, n_process=125)
    # Set param values
    α = 2.5e-1

    # Simulate true power-law adaptation
    @info "Evaluating real power-law"
    @time y_pla, I_pla = adapt_pla(x, α, 5e-3; fs=100e3)

    # Simulate parallel set of exponential processes
    @info "Evaluating iir power-law approximation"
    @time y_ea, I_ea, I_mat = adapt_ea_iir_parallel(x, 0.1 * α, LogRange(τ_short, τ_long, n_process); fs=100e3)

    # Plot
    fig = Figure(; resolution=(900, 700))
    axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
    ax_time = Axis(fig[3, 1])
    ax_compare = Axis(fig[4, :])
    lines!(axs[1, 1], t, I_pla)
    lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
    lines!(axs[2, 1], t, I_ea)
    lines!(axs[1, 2], t, y_pla)
    lines!(axs[2, 2], t, y_ea)
    cols = collect(eachcol(I_mat))[1:5:end]
    colors = get(ColorSchemes.batlow, LinRange(0.0, 1.0, length(cols)))
    map(zip(colors, cols)) do (color, col)
        lines!(ax_time, t, col; color=color)
    end
    ylims!.(axs[:, 1], 0.0, 4000.0)
    lines!(ax_compare, y_pla .- y_ea)
    fig

    # Check tolerance
    rtol = 0.0
    match = false
    while !match
        rtol += 0.01
        if isapprox(y_pla, y_ea; rtol=rtol)
            match = true
        end
    end
    @info "Approximate match at rtol = $rtol"
end
fig2()
