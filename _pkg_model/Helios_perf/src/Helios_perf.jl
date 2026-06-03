module Helios_perf

export format_axis_time_time!
export c_baseline
export p_performance, p_performance_prealloc, p_performance_prealloc_benefit, p_channel_scaling, p_performance_timevec

using BenchmarkTools
using CairoMakie
using Dates
using Helios

# ##########################################################################################
# Define functions to format axes as needed
# ##########################################################################################
function format_axis_time_time!(ax; t_min=1e-2, t_max=1e2)
    xlims!(ax, t_min, t_max)
    ylims!(ax, t_min, t_max)
    ax.xscale = log10
    ax.yscale = log10
    t_min = 10^floor(log10(t_min))
    t_max = 10^ceil(log10(t_max))
    n_tick = Int(log10(t_max) - log10(t_min) + 1)
    ax.xticks = 10 .^ (LinRange(log10(t_min), log10(t_max), n_tick))
    ax.yticks = 10 .^ (LinRange(log10(t_min), log10(t_max), n_tick))
    ax.xminorticksvisible = true
    ax.xminorticks = IntervalsBetween(9)
    ax.yminorticksvisible = true
    ax.yminorticks = IntervalsBetween(9)
    ax.xlabel = "Simulation time (s)"
    ax.ylabel = "Execution time (s)"
    return ax
end

function format_axis_channels_time!(ax; t_min=1e-2, t_max=1e2, chan_min=1, chan_max=128)
    xlims!(ax, chan_min/2, chan_max*2)
    ylims!(ax, t_min, t_max)
    ax.xscale = log10
    ax.yscale = log10
    t_min = 10^floor(log10(t_min))
    t_max = 10^ceil(log10(t_max))
    n_tick = Int(log10(t_max) - log10(t_min) + 1)
    ax.yticks = 10 .^ (LinRange(log10(t_min), log10(t_max), n_tick))
    chan_min = 2^floor(log2(chan_min))
    chan_max = 2^ceil(log2(chan_max))
    n_tick = Int(log2(chan_max) - log2(chan_min) + 1)
    ax.xticks = 2 .^ (LinRange(log2(chan_min), log2(chan_max), n_tick))
    ax.yminorticksvisible = true
    ax.yminorticks = IntervalsBetween(9)
    ax.xlabel = "Number of channels"
    ax.ylabel = "Execution time (s)"
    return ax
end

# ##########################################################################################
# Define functions to generate dummy inputs, prefixed with d_
# ##########################################################################################
function d_stim(n=10000)
    zeros(n)
end

function d_cf(n=1)
    if n == 1
        return [1e3]
    else
        exp.(LinRange(log(0.125e3), log(16e3), n))
    end
end

# ##########################################################################################
# Define a number of functions that implement example calls to Helios functions for testing
# Each function is prefixed with c_
# ##########################################################################################
function c_baseline()
    sim_gfc2023(d_stim(), d_cf());
end

function c_multichan_small()
    sim_gfc2023(d_stim(), d_cf(5));
end

function c_freeze_gain()
    sim_gfc2023(d_stim(), d_cf(); moc_fix_gain=true)
end

function c_time(n)
    sim_gfc2023(d_stim(n), d_cf());
end

function c_prealloc(len_total=20000, n_chan=101)
    GFC2023_Mem(len_total, n_chan)
end

# ##########################################################################################
# Define functions to generate automatic plots of performance data
# Each function is prefixed with p_
# ##########################################################################################
const p_path = "C:\\Users\\dguest2\\cl_code\\Helios\\Helios_perf\\figures"

# Date-tagged plot saving function
function p_save(fig, fn)
    folder = joinpath(p_path, string(today()))
    if !isdir(folder)
        mkdir(folder)
    end
    save(joinpath(folder, fn), fig)
end

# Create scatterviolins of performance data for different reference c_funcs
function p_performance()
    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Function name", ylabel="Time (s)", yscale=log10)

    # Select which c_funcs to test
    c_funcs = [c_baseline, c_multichan_small, c_freeze_gain]

    # Loop through c_funcs and run each 100 times
    map(enumerate(c_funcs)) do (idx, c_func)
        # Print progress
        println("Testing $(c_func)")

        # Measure execution time
        times = map(1:100) do _ 
            @elapsed c_func()
        end

        # Scatter the raw time recordings
        scatter!(ax, fill(idx, length(times)), times, markersize=5, color=:black)

        # Mark the minimum and median times
        min_time = minimum(times)
        median_time = median(times)
        scatter!(ax, [idx+0.1], [min_time], color=:black, markersize=10.0, marker='<')
        scatter!(ax, [idx+0.1], [median_time], color=:black, markersize=10.0, marker='<')
    end

    # Adjust xticks to indicate cfunc values
    ax.xticks = (1:length(c_funcs), string.(c_funcs))

    # Adjust yticks to be appropriate for yscale
    ax.yticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]
    ax.yminorticksvisible = true
    ax.yminorticks = IntervalsBetween(9)
    ylims!(ax, 1e-3, 1e2)
    xlims!(ax, 0.5, length(c_funcs) + 0.5)

    # Save figure
    p_save(fig, "01_overall_performance.png")
end

# Create function to plot performance vs time vector length
function p_performance_timevec(n_rep=25)
    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1])
    format_axis_time_time!(ax; t_min=1e-3, t_max=1e1)

    # Select time vector lengths to test    
    Ns = 1000 .* 2 .^ (1:8)

    # Loop through c_funcs and run each 100 times
    mins = map(enumerate(Ns)) do (idx, n)
        # Print progress
        println("Testing n = $(n)")

        # Measure execution time
        times = map(1:n_rep) do _ 
            @elapsed c_time(n)
        end

        # Scatter the raw time recordings
        scatter!(ax, fill(n/100e3, length(times)), times, markersize=5, color=:black)

        # Mark the minimum and median times
        min_time = minimum(times)
        median_time = median(times)
        scatter!(ax, [n/100e3*1.20], [min_time], color=:black, markersize=10.0, marker='<')
        scatter!(ax, [n/100e3*1.20], [median_time], color=:black, markersize=10.0, marker='<')
        return min_time
    end

    # Add predicted linear scaling from smallest # of samples 
    # lines!(ax, Ns ./ 100e3, mins[1] .* (Ns ./ Ns[1]) ./ 100e3; linestyle=:dash, color=:red)

    # Add predicted linear scaling from 16000 samples (assuming 16000 is 4)
    # lines!(ax, Ns ./ 100e3, mins[4] .* (Ns ./ Ns[4]) ./ 100e3; linestyle=:dash, color=:blue)
    x = exp.(LinRange(log(1e-3), log(1e1), 1000))
    y = x
    lines!(ax, x, y; color=:black, linestyle=:dash)
    lines!(ax, x, y .* 10; color=:gray, linestyle=:dash)
    lines!(ax, x, y ./ 10; color=:gray, linestyle=:dash)

    # Save figure
    p_save(fig, "06_scaling_time.png")
    fig
end

# Create scatterviolins of performance data for different reference c_funcs focusing on preallocation issues
function p_performance_prealloc(n_rep=20; fs=100e3)
    # Create figure
    fig = Figure(; size=(1200, 400))

    # Loop through c_funcs and run each 100 times
    Ns = [2500, 5000, 10000, 20000, 40000, 80000]
    n_chans = [1, 5, 21, 101]
    for (idx, n_chan) in enumerate(n_chans)
        # Create axis 
        ax = Axis(fig[1, idx])
        format_axis_time_time!(ax; t_min=1e-3, t_max=1e1)
        ax.title = "n_chan = $(n_chan)"
        ax.xlabel = "Allocated memory duration (s)"
        if idx > 1
            ax.ylabel = ""
            ax.yticklabelsvisible = false
        end

        # Loop over lengths
        for n in Ns
            # Print progress
            println("Testing $(n) samples...")

            # Measure execution time
            times = map(1:n_rep) do _ 
                @elapsed c_prealloc(n, n_chan)
            end

            # Scatter the raw time recordings
            scatter!(ax, fill(n/fs, length(times)), times, markersize=5, color=:black)

            # Mark the minimum and median times
            min_time = minimum(times)
            median_time = median(times)
            scatter!(ax, [n/fs*1.1], [min_time], color=:black, markersize=10.0, marker='<')
            scatter!(ax, [n/fs*1.1], [median_time], color=:black, markersize=10.0, marker='<')
        end
    end

    # Save figure
    p_save(fig, "02_allocation_performance.png")
    fig
end

function p_performance_prealloc_benefit(n_rep=25; t_min=1e-2, t_max=1e1)
    # Select range of number of channels and repetitions
    N_chans = [1, 4, 16, 64]
    durs = exp.(LinRange(log(5e-2), log(1e0), 9))

    # Create figure
    fig = Figure(; size=(400*length(N_chans), 400))

    # Loop over number of channels
    for (idx, n_chan) in enumerate(N_chans)
        # Create axis
        ax = Axis(fig[1, idx])
        ax = format_axis_time_time!(ax; t_min=t_min, t_max=t_max)
        ax.title = "n_chan = $(n_chan)"

        # Loop over durations
        for dur in durs
            # Print progress
            println("Testing n_chan = $(n_chan)... dur = $(dur)...")

            # Set parameters
            y = d_stim(Int(round(dur * 100e3)))
            cfs = d_cf(n_chan)

            # Set up pre-allocated memory
            mem = GFC2023_Mem(y, cfs; dur_pad_left=0.0, dur_pad_right=0.0)

            # Measure regular sim_gfc2023 calls several times
            t_reg = map(1:n_rep) do _
                @elapsed sim_gfc2023(y, cfs; fractional=false)
            end
            scatter!(ax, fill(dur, length(t_reg)), t_reg; color=:black)
            scatter!(ax, dur*1.1, minimum(t_reg); color=:black, markersize=10.0, marker='<')
            scatter!(ax, dur*1.1, median(t_reg); color=:black, markersize=10.0, marker='<')

            # Measure sim_gfc2023 calls with pre-allocated memory several times
            t_prealloc = map(1:n_rep) do _
                @elapsed sim_gfc2023!(mem, y, cfs; dur_pad_left=0.0, dur_pad_right=0.0)
            end
            scatter!(ax, fill(dur, length(t_prealloc)), t_prealloc; color=:red)
            scatter!(ax, dur*1.1, minimum(t_prealloc); color=:red, markersize=10.0, marker='<')
            scatter!(ax, dur*1.1, median(t_prealloc); color=:red, markersize=10.0, marker='<')
        end
    end

    # Return fig
    fig
end

function p_channel_scaling(n_rep=40; moc_width=0.8)
    # Select numbers of channels to test
    N = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1])
    format_axis_channels_time!(ax; t_min=1e-3, t_max=1e2, chan_min=1, chan_max=maximum(N))

    # Test different numbers of channels
    meds = map(enumerate(N)) do (idx, n)
        # Print progress
        println("Running n = $(n) channels...")

        # Measure execution time for different numbers of channels
        times = map(1:n_rep) do _
            @elapsed sim_gfc2023(d_stim(), d_cf(n); moc_width=0.8)
        end

        # Plot the results
        scatter!(ax, fill(float(n), length(times)), times; color=:black)

        # Return median
        return median(times)
    end

    # Add scaling from t_min*10 as linear function of N
    x = 2 .^ LinRange(log2(1.0), log2(maximum(N)), 1000)
    y = x .* 10 .* 1e-3
    lines!(ax, x, y; color=:black, linestyle=:dash)
    lines!(ax, x, y .* 10; color=:gray, linestyle=:dash)
    lines!(ax, x, y ./ 10; color=:gray, linestyle=:dash)

    # Save figure and return 
    p_save(fig, "03_channel_scaling_moc_width=$(moc_width).png")
    return fig
end

end # module Helios_perf
