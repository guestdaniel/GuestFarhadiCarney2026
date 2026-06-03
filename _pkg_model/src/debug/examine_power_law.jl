using Helios
using CairoMakie
using AuditorySignalUtils
using DSP
using Optim

# Simulate response to 1-kHz pure tone at 50 dB SPL
stim = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 1.0, 100e3), 0.01, 100e3), 50.0)
resp = sim_gfc2023_dict(stim, 1000.0; moc_weight_ic=0.0, moc_weight_wdr=0.0, powerlaw_mode=1)
resp_old = sim_orig_dict(stim, 1000.0)
x = resample(resp["expon"], 10/100)
x = vcat(x, zeros(5000))
t = 0.0:(1/10e3):(length(x)/10e3-1/10e3)

# Function to implement PLA in Julia
function pla(x, α, β)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * 1/10e3 / ((n-j)*1/10e3 + β)
        end
    end
    return y, I
end

function ea(x, τₐ, τₑ)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * 1/10e3 * exp( ((j-n) * (1/10e3)) * (1/τₑ))
        end
    end
    return y, I
end

function ea_hack(x, τₐ, τ)
    y = zeros(size(x))
    I = zeros(size(x))
    d = exp(-(1/10e3)/τ)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        if n == 1
            I[n] = (1-d) * 1e-1 * y[n]
        else
            I[n] = (1-d) * 1e-1 * y[n] + d * I[n-1]
        end
    end
    return y, I
end

function ea_parallel(x, τₐ, τs)
    y = zeros(length(x))
    I = zeros(length(x), length(τs))
    I_comb = zeros(length(x))
    ds = @. exp(-(1/10e3)/τs)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I_comb[n-1])
        end

        # Compute I[n]
        for j in eachindex(ds)
            if n == 1
                I[n, j] = (1-ds[j]) * 1e-1 * y[n]
            else
                I[n, j] = (1-ds[j]) * 1e-1 * y[n] + ds[j] * I[n-1, j]
            end
        end
        I_comb[n] = sum(I[n, :])
    end
    return y, I_comb
end

## Figure #1: Exponential adaptation via direct integration vs IIR EWMA 

# Simulate direct integration
y_direct, I_direct = ea(x, 1e-1, 1e-1)
# Simulate IIR EWMA
y_iir, I_iir = ea_hack(x, 1e-1, 1e-1)

# Plot
fig = Figure(; resolution=(900, 500))
axs = [Axis(fig[i, j]) for i in 1:2, j in 1:3]
axs_compare = [Axis(fig[3, 2]), Axis(fig[3, 3])]
lines!(axs[1, 1], t, x)
lines!(axs[1, 2], t, I_direct)
lines!(axs[1, 3], t, y_direct)
lines!(axs[2, 1], t, x)
lines!(axs[2, 2], t, I_iir)
lines!(axs[2, 3], t, y_iir)
lines!(axs_compare[1], t, I_direct)
lines!(axs_compare[1], t, I_iir)
lines!(axs_compare[2], t, y_direct)
lines!(axs_compare[2], t, y_iir)
fig

## Figure #2: Compare adaptation integrals and outputs for IIR parallel vs real PLA

# Simulate true power-law adaptation
@time y_pla, I_pla = pla(x, 2.5e-6*100e3, 5e-3)

# Simulate parallel set of exponential processes
@time y_ea, I_ea = ea_parallel(x, 0.7e-1, LogRange(5e-3, 1e4, 100))
@time y_ea, I_ea = ea_parallel(x, 25 * 2.5e-6*100e3, LogRange(5e-3, 1e4, 200))

# Plot
fig = Figure(; resolution=(900, 400))
axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
ax_compare = Axis(fig[3, :])
lines!(axs[1, 1], t, I_pla)
lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
lines!(axs[2, 1], t, I_ea)
lines!(axs[1, 2], t, y_pla)
lines!(axs[2, 2], t, y_ea)
ylims!.(axs[:, 1], 0.0, 4000.0)
lines!(ax_compare, y_pla .- y_ea)
fig

## Figure #3: Optimize scale factor, upper/lower time constants, and number of processes
# Select range option uphere
rangefunc = LinRange
# Define our functions
function foptim(x) 
    error = ea_parallel(x, x[1], rangefunc(exp(x[2]), exp(x[3]), Int(round(max(0.0, x[4])))))[1] .- pla(x, 2.5e-6*100e3, 5e-3)[1]
    return rms(error)
end

function fdisplay(x̂; xlims=(0.0, 0.1)) 
    y_ea, I_ea = ea_parallel(x, x̂[1], rangefunc(exp(x̂[2]), exp(x̂[3]), Int(round(x̂[4]))))
    y_pla, I_pla = pla(x, 2.5e-6*100e3, 5e-3)
    fig = Figure(; resolution=(1500, 900))
    axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
    ax_compare = Axis(fig[3, :])
    lines!(axs[1, 1], t, I_pla)
    lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
    lines!(axs[2, 1], t, I_ea)
    lines!(axs[1, 2], t, y_pla)
    lines!(axs[2, 2], t, y_ea)
    ylims!.(axs[:, 1], 0.0, 4000.0)
    lines!(ax_compare, t, y_pla .- y_ea)
    xlims!.(axs, xlims...)
    xlims!(ax_compare, xlims...)
    ax_compare.title = "RMS error = $(rms(y_pla .- y_ea))"
    fig
end

function fprint(x)
    xtemp = copy(x)
    xtemp[2:3] = exp.(xtemp[2:3])
    xtemp = round.(xtemp; digits=5)
    "Optimized values are: Coefficient of $(xtemp[1]), time-constant range of ($(xtemp[2]), $(xtemp[3])) s, process count of $(xtemp[4])"
end

x_init = [25*2.5e-6*100e3, log(5e-3), log(1e4), 100]

results = optimize(foptim, x_init)
x̂ = Optim.minimizer(results)
fprint(x̂)
fdisplay(x̂; xlims=(0.0, 1.0))

## Figure #4: Optimize scale factor, upper/lower time constants, and number of processes
## Use arbitrary scaling for range
# Select range option uphere
customrange(α) = (x...) -> LinRange(x[1] ^ (1/α), x[2] ^ (1/α), x[3]) .^ α
# Define our functions
function foptim(x) 
    y_ea = ea_parallel(x, x[1], customrange(exp(x[5]))(exp(x[2]), exp(x[3]), Int(round(max(0.0, x[4])))))[1] 
    y_pla = pla(x, 2.5e-6*100e3, 5e-3)[1]
    return rms(y_pla .- y_ea)
end

function fdisplay(x̂; xlims=(0.0, 0.1)) 
    y_ea, I_ea = ea_parallel(x, x̂[1], customrange(exp(x[5]))(exp(x̂[2]), exp(x̂[3]), Int(round(x̂[4]))))
    y_pla, I_pla = pla(x, 2.5e-6*100e3, 5e-3)
    fig = Figure(; resolution=(1500, 900))
    axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
    ax_compare = Axis(fig[3, :])
    lines!(axs[1, 1], t, I_pla)
    lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
    lines!(axs[2, 1], t, I_ea)
    lines!(axs[1, 2], t, y_pla)
    lines!(axs[2, 2], t, y_ea)
    ylims!.(axs[:, 1], 0.0, 4000.0)
    lines!(ax_compare, t, y_pla .- y_ea)
    xlims!.(axs, xlims...)
    xlims!(ax_compare, xlims...)
    ax_compare.title = "RMS error = $(rms(y_pla .- y_ea))"
    fig
end

function fprint(x)
    xtemp = copy(x)
    xtemp[2:3] = exp.(xtemp[2:3])
    xtemp[5] = exp(xtemp[5])
    xtemp = round.(xtemp; digits=5)
    "Optimized values are: Coefficient of $(xtemp[1]), time-constant range of ($(xtemp[2]), $(xtemp[3])) s, process count of $(xtemp[4]), range exponent of $(xtemp[5])"
end

x_init = [25*2.5e-6*100e3, log(5e-3), log(1e4), 100, log(1.0)]

results = optimize(foptim, x_init)
x̂ = Optim.minimizer(results)
fprint(x̂)
fdisplay(x̂; xlims=(0.0, 1.0))

## Figure #4: Optimize scale factor, upper/lower time constants, and number of processes
## Use arbitrary scaling for range
# Select range option uphere
customrange(α) = (x...) -> LinRange(x[1] ^ (1/α), x[2] ^ (1/α), x[3]) .^ α
lp = digitalfilter(Lowpass(200.0; fs=10e3), Butterworth(4))
lpf(x) = filt(lp, x)

# Define our functions
function foptim(x) 
    y_ea = ea_parallel(x, exp(x[1]), customrange(exp(x[5]))(exp(x[2]), exp(x[3]), Int(round(max(0.0, x[4])))))[1]
    y_pla = pla(x, 2.5e-6*100e3, 5e-3)[1]
    return rms(lpf(y_pla) .- lpf(y_ea))
end

function fdisplay(x̂; xlims=(0.0, 0.1)) 
    y_ea, I_ea = ea_parallel(x, exp(x̂[1]), customrange(exp(x[5]))(exp(x̂[2]), exp(x̂[3]), Int(round(x̂[4]))))
    y_pla, I_pla = pla(x, 2.5e-6*100e3, 5e-3)
    fig = Figure(; resolution=(1500, 900))
    axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
    ax_compare = Axis(fig[3, :])
    lines!(axs[1, 1], t, I_pla)
    lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
    lines!(axs[2, 1], t, I_ea)
    lines!(axs[1, 2], t, y_pla)
    lines!(axs[2, 2], t, y_ea)
    ylims!.(axs[:, 1], 0.0, 4000.0)
    lines!(ax_compare, t, lpf(y_pla) .- lpf(y_ea))
    xlims!.(axs, xlims...)
    xlims!(ax_compare, xlims...)
    ax_compare.title = "RMS error = $(rms(lpf(y_pla) .- lpf(y_ea)))"
    fig
end

function fprint(x)
    xtemp = copy(x)
    xtemp[1] = exp(xtemp[1])
    xtemp[2:3] = exp.(xtemp[2:3])
    xtemp[5] = exp(xtemp[5])
    xtemp = round.(xtemp; digits=5)
    "Optimized values are: Coefficient of $(xtemp[1]), time-constant range of ($(xtemp[2]), $(xtemp[3])) s, process count of $(xtemp[4]), range exponent of $(xtemp[5])"
end

x_init = [log(25*2.5e-6*100e3), log(5e-3), log(1e4), 100, log(2.0)]
x_low = [log(1e-3), log(1e-5), log(1e1), 10, log(1.0)]
x_high = [log(1e3), log(1e0), log(1e10), 1000, log(1000.0)]

results = optimize(foptim, x_low, x_high, x_init, SAMIN())
x̂ = Optim.minimizer(results)
fprint(x̂)
fdisplay(x̂; xlims=(0.0, 1.5))

## Figure #5: Back to logscale
# Select range option uphere
lp = digitalfilter(Lowpass(200.0; fs=10e3), Butterworth(4))
lpf(x) = filt(lp, x)

# Define our functions
function foptim(x) 
    y_ea = ea_parallel(x, exp(x[1]), LogRange(exp(x[2]), exp(x[3]), Int(round(max(0.0, x[4])))))[1]
    y_pla = pla(x, 2.5e-6*100e3, 5e-3)[1]
    return rms(lpf(y_pla) .- lpf(y_ea))
end

function fdisplay(x̂; xlims=(0.0, 0.1)) 
    y_ea, I_ea = ea_parallel(x, exp(x̂[1]), LogRange(exp(x̂[2]), exp(x̂[3]), Int(round(x̂[4]))))
    y_pla, I_pla = pla(x, 2.5e-6*100e3, 5e-3)
    fig = Figure(; resolution=(1500, 900))
    axs = [Axis(fig[i, j]) for i in 1:2, j in 1:2]
    ax_compare = Axis(fig[3, :])
    lines!(axs[1, 1], t, I_pla)
    lines!(axs[2, 1], t, I_pla; color=:pink, linewidth=0.5)
    lines!(axs[2, 1], t, I_ea)
    lines!(axs[1, 2], t, y_pla)
    lines!(axs[2, 2], t, y_ea)
    ylims!.(axs[:, 1], 0.0, 4000.0)
    lines!(ax_compare, t, lpf(y_pla) .- lpf(y_ea))
    xlims!.(axs, xlims...)
    xlims!(ax_compare, xlims...)
    ax_compare.title = "RMS error = $(rms(lpf(y_pla) .- lpf(y_ea)))"
    fig
end

function fprint(x)
    xtemp = copy(x)
    xtemp[1] = exp(xtemp[1])
    xtemp[2:3] = exp.(xtemp[2:3])
    xtemp = round.(xtemp; digits=5)
    "Optimized values are: Coefficient of $(xtemp[1]), time-constant range of ($(xtemp[2]), $(xtemp[3])) s, process count of $(xtemp[4])"
end

x_init = [log(25*2.5e-6*100e3), log(5e-3), log(1e4), 100]
x_low = [log(1e-3), log(1e-5), log(1e1), 10]
x_high = [log(1e3), log(1e0), log(1e10), 1000]

results = optimize(foptim, x_low, x_high, x_init, SAMIN())
x̂ = Optim.minimizer(results)
fprint(x̂)
fdisplay(x̂; xlims=(0.0, 1.5))

