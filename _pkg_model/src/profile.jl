export measure_linear_channel_scaling, measure_current_vs_2014, measure_all

"""
    test_single_channel(n)

Tests n reps of a single channel response using `sim_gfc2023`
"""
function test_single_channel(n_rep; dur=1.0, fs=100e3, kwargs...)
    stim = zeros(Int(round(dur*fs)))
    for _ = 1:n_rep
        sim_gfc2023(stim, 1000.0; powerlaw_mode=0, fs=fs, kwargs...)
    end
end

"""
    test_multi_channel(n)

Tests n reps of a multi channel response using `sim_gfc2023`
"""
function test_multi_channel(n_chan, n_rep; dur=1.0, fs=100e3, kwargs...)
    stim = zeros(Int(round(dur*fs)))
    cf = LogRange(0.2e3, 20e3, n_chan)
    for _ = 1:n_rep
        sim_gfc2023(stim, cf; powerlaw_mode=0, fs=fs, kwargs...)
    end
end

"""
    test_single_channel_2014(n)

Tests n reps of a single channel response using `sim_anrate_zbc2014`
"""
function test_single_channel_2014(n_rep; dur=1.0, fs=100e3)
    stim = zeros(Int(round(dur*fs)))
    for _ = 1:n_rep
        ihc = sim_ihc_zbc2014(stim, 1000.0)
        sim_anrate_zbc2014(ihc, 1000.0)
    end
end

"""
    profile_single_channel()

Profile simulations of several 1-s single-channel simulation runs for the efferent model
"""
function profile_single_channel(; n_rep=10, dur=1.0)
    Profile.clear()
    @profile test_single_channel(n_rep; dur=dur)
    ProfileView.view(; C=true)
end

"""
    profile_single_channel_2014()

Profile simulations of several 1-s single-channel simulation runs for the 2014 model
"""
function profile_single_channel_2014(; n_rep=10, dur=1.0)
    Profile.clear()
    @profile test_single_channel_2014(n_rep; dur=dur)
    ProfileView.view(; C=true)
end

"""
    profile_multi_channel()

Profile simulations of several 1-s multi-channel simulation runs for the efferent model using
"""
function profile_multi_channel(; n_chan=10, n_rep=10, dur=1.0)
    Profile.clear()
    @profile test_multi_channel(n_chan, n_rep; dur=dur)
    ProfileView.view(; C=true)
end

"""
    measure_linear_channel_scaling()

We want to ensure that adding additional channels simply results in a proportional increase 
in run time (i.e., there is no fundamental overhead associated with running multi-channel
simulations). Here, we compare theoretical predictions of multi-channel runtime based on 
single-channel run time to empirical multi-channel runtime.
"""
function measure_linear_channel_scaling(
    n_chans=(2 .^ (0:1:5));
    n_rep=3,
)
    # Calculate runtime for single run, and then extrapolate multi-channel runs
    runtime_single = (@elapsed test_single_channel(n_rep))/n_rep
    runtime_single = n_chans .* runtime_single

    # Empirically measure single-channel runtime
    runtime_single_em = map(n_chans) do n_chan
        (@elapsed test_single_channel(n_rep*n_chan))/n_rep
    end

    # Empirically measure multi-channel runtimes
    runtime_multi = map(n_chans) do n_chan
        (@elapsed test_multi_channel(n_chan, n_rep))/n_rep
    end

    # Empirically measure multi-channel runtimes w/o cross-frequency connections
    runtime_multi_no_cross = map(n_chans) do n_chan
        (@elapsed test_multi_channel(n_chan, n_rep; moc_width_wdr=0.0))/n_rep
    end

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, yminorticksvisible=true)

    # Plot and label each set of data
    lines!(ax, n_chans, runtime_single; color=:gray, linestyle=:dash, label="Single channel (extrapolated)")
    scatter!(ax, n_chans, runtime_single; color=:gray)

    lines!(ax, n_chans, runtime_single_em; color=:black, label="Single channel (measured)")
    scatter!(ax, n_chans, runtime_single_em; color=:black)

    lines!(ax, n_chans, runtime_multi; color=:red, label="Multi channel")
    scatter!(ax, n_chans, runtime_multi; color=:red)

    lines!(ax, n_chans, runtime_multi_no_cross; color=:blue, label="Multi channel (width=0.0)")
    scatter!(ax, n_chans, runtime_multi; color=:blue)

    # Add label denoting performance gap at largest tests number of channels
    lines!(ax, 1.3 .* [n_chans[end], n_chans[end]], [runtime_single_em[end], runtime_multi_no_cross[end]]; color=:purple)
    gap = runtime_multi_no_cross[end]/runtime_single_em[end]
    text!(ax, [1.35 * n_chans[end]], [exp((1/2) * (log(runtime_single_em[end]) + log(runtime_multi_no_cross[end])))]; text="$(round(gap; digits=2))", color=:purple)

    # Label figure
    ax.xlabel = "Number of channels"
    ax.ylabel = "Average runtime (s)"
    ax.title = "Linear channel scaling; n_rep=$n_rep"
    ax.yticks = [0.1, 1.0, 10.0]
    ax.xticks = n_chans
    ax.yminorticks = IntervalsBetween(9)
    axislegend(ax; position=:rb)
    xlims!(ax, 0.5 * n_chans[1], 2.0 * n_chans[end])
    fig

    # Auto-save figure with date tag
    save(projectdir("plots", "$(today())_measure_linear_channel_scaling.png"), fig)

    # Return figure
    fig
end

"""

    measure_current_vs_2014()

It is important to know how runtimes in the original 2014 model stack up against runtimes 
in the efferent model under various conditions. This function measures single and 
multi-channel simulation performance for the original 2014 model code (or, at least, my 
adapted Julia version) versus the current efferent model code. 
"""
function measure_current_vs_2014(
    durs=LogRange(0.1, 10.0, 9);
    n_rep=10,
)
    # Measure runtime for 2014 model
    runtimes_2014 = map(durs) do dur
        (@elapsed test_single_channel_2014(n_rep; dur=dur))/n_rep
    end

    runtimes = map(durs) do dur
        (@elapsed test_single_channel(n_rep; dur=dur))/n_rep
    end

    # Plot results
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, yminorticksvisible=true, xminorticksvisible=true)
    lines!(ax, durs, runtimes_2014; color=:gray, label="2014")
    scatter!(ax, durs, runtimes_2014; color=:gray)
    lines!(ax, durs, runtimes; color=:red, label="Efferent")
    scatter!(ax, durs, runtimes; color=:red)
    ax.xlabel = "Simulation duration (s)"
    ax.ylabel = "Average runtime (s)"
    ax.title = "Efferent model performance; n_rep=$n_rep"
    ax.yticks = 1.0 .* 10.0 .^ (-4:1:2)
    ax.xticks = 1.0 .* 10.0 .^ (-2:1:2)
    xlims!(ax, 1e-2, 1e2)
    ylims!(ax, 1e-4, 1e2)
    ax.yminorticks = IntervalsBetween(9)
    ax.xminorticks = IntervalsBetween(9)
    axislegend(ax)
    fig

    # Auto-save figure with date tag
    save(projectdir("plots", "$(today())_measure_current_vs_2014.png"), fig)

    # Return figure
    fig
end

"""
    profile_single_vs_multi()

Do a profiling run for total-sample-matched single vs multichannel runs
"""
function profile_single_vs_multi()
    profile_single_channel(; n_rep=10)
    profile_multi_channel(; n_rep=1, n_chan=10)
end

"""
    measure_all()

Runs all measurement functions to generate a fresh set of figures 
"""
function measure_all()
    measure_linear_channel_scaling()
    measure_current_vs_2014()
end
