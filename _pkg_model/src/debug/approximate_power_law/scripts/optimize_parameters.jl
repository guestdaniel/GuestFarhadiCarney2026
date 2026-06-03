using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using ColorSchemes
using Optim
using Printf

## Define function to assist with optimization 
function adaptation_factory(x; α=2.5e-1, β=5e-4, fs_orig=100e3, fs=100e3, n_process=125)
    # Zero-pad
    x = vcat(x, zeros(Int(round(0.1*fs_orig))))

    # Compute model response
    resp = sim_gfc2023_dict(x, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, fs=fs_orig)

    # Extract and downsample output waveform from exponential adaptation stage (input to PLA)
    expon = resample(resp["expon"], fs/fs_orig)

    # Simulate true power-law adaptation
    y_pla = adapt_pla_clike(expon, α, β; fs=fs)
    y_internal = resample(resp["sout1"], fs/fs_orig)

    # Create time-axis
    t = 0.0:(1/fs):(length(y_pla)/fs - 1/fs)

    # Emit functions for optimization using 125 parallel processes
    function f(θ)
        y_ea, _, _ = adapt_ea_iir_parallel(expon, θ[1] * α, LogRange(exp(θ[2]), exp(θ[3]), n_process); fs=fs)
        rms(y_pla .- y_ea)
    end

    return t, expon, y_pla, y_internal, f
end

## Optimize parameters in response to 1-kHz pure tone
stim = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.2, 100e3), 0.01, 100e3), 50.0)
α = 2.5e-1
β = 5e-4
n_process = 200
x_init = [1.0, log(5e-4), log(1e4)]
t, expon, y_pla, y_internal, f = adaptation_factory(stim; α=α, β=β, fs=100e3, n_process=n_process)
results = optimize(f, x_init)
x̂ = Optim.minimizer(results)
ŷ, _, _ = adapt_ea_iir_parallel(expon, x̂[1] * alpha, LogRange(exp(x̂[2]), exp(x̂[3]), n_process); fs=100e3)

## Plot results of optimization
fig = Figure(; resolution=(1300, 1000))

# Plot real vs approximate PLA
ax = Axis(fig[1, 1])
ax.title = "True PLA (black) vs approximate PLA (pink)"
lines!(ax, t, y_pla; color=:black)
lines!(ax, t, ŷ; color=:pink)
xlims!(0.0, 0.3)
ax = Axis(fig[1, 2])
ax.title = "True PLA − approximate PLA"
lines!(ax, t, y_internal .- ŷ; color=:black)
xlims!(0.0, 0.3)

# Plot real vs approximate PLA
ax = Axis(fig[2, 1])
ax.title = "True PLA (black) vs approximate PLA (pink)"
lines!(ax, t, y_pla; color=:black)
lines!(ax, t, ŷ; color=:pink)
xlims!(0.0, 0.01)
ax = Axis(fig[2, 2])
ax.title = "True PLA − approximate PLA"
lines!(ax, t, y_internal .- ŷ; color=:black)
xlims!(0.0, 0.01)

# Label
Label(fig[0, :], "$n_process processes, RMS error = $(rms(y_internal .- ŷ))")
fig