using Distributed
addprocs(6)
@everywhere using Pkg
@everywhere Pkg.activate("TestBench")
@everywhere begin
    using Helios
    using Statistics
end
using CairoMakie
using UtilitiesViz
using Colors
using ColorSchemes

## Set parameters
function test(
    dur=0.5,
    fs=100e3,
    n_rep=10,
    stage="hsr",
)
    # Set stimulus
    stim = zeros(Int(round(dur*fs)))

    # Map over interval averaging vs no averaging
    results_no_vs_internal = map([sim_gfc2023_dict, sim_gfc2023ia_dict]) do func
        # Map over n_rep
        pmap(1:n_rep) do n
            println("Doing rep $n of $n_rep for $(string(func))")
            resp = func(stim, [1000.0]; fractional=true, powerlaw_mode=2)
            resp[stage][1]
        end
    end

    # Get external averaging results
    external = pmap(1:n_rep) do n
        println("Doing rep $n of $n_rep for external averaging")
        n_sim_avg = stage == "lsr" ? 5 : 10
        resps = map(1:n_sim_avg) do _
            resp = sim_gfc2023_dict(stim, [1000.0]; fractional=true, powerlaw_mode=2)
            resp[stage][1]
        end
        mean(resps)
    end

    return results_no_vs_internal, external
end

function doplot_waveform_examples(stage="hsr", ylims=([0.0, 350.0], [0.0, 75.0]))
    # Collect data
    (noavg, internal), external = test(0.5, 100e3, 10, stage)
    results = [noavg, external, internal]

    ## Plot
    # Create variables we need
    t = 0.0:(1/100e3):(length(results[1][1])/100e3-1/100e3)
    results = [noavg, external, internal]
    titles = ["No averaging", "External averaging", "Internal averaging"]
    colors = ColorSchemes.tab10

    # Create figures and axs
    fig = Figure(; resolution=(1100, 1500))
    axs_time = [Axis(fig[i, 1]) for i in eachindex(results)]
    axs_hist = [Axis(fig[i, 2]) for i in eachindex(results)]
    axs_μ = [Axis(fig[i, 3]) for i in eachindex(results)]
    axs_σ = [Axis(fig[i, 4]) for i in eachindex(results)]
    axs = hcat(axs_time, axs_hist, axs_μ, axs_σ)

    # Plot time-domain waveforms
    map(enumerate(zip(axs_time, titles, results))) do (idx, (ax, title, result))
        map(zip(result, colors)) do (r, c)
            lines!(ax, t, r; color=c)
        end
        ax.title = title
        ax.xlabel = "Time (s)"
    end

    # Plot density kernel histograms of instantaneous rates
    map(enumerate(zip(axs_hist, titles, results))) do (idx, (ax, title, result))
        map(zip(result, colors)) do (r, c)
            density!(ax, r[t .> 0.1]; direction=:y, color=(c, 0.3), strokecolor=c, strokewidth=2.0, strokearound=true)
        end
    end

    # Plot summary stats for mean rates
    map(enumerate(zip(axs_μ, titles, results))) do (idx, (ax, title, result))
        M = map(zip(result, colors)) do (r, c)
            μ = mean(r[t .> 0.1])
            scatter!(ax, [1.0], [μ]; color=c)
            μ
        end
        scatter!(ax, [2.0], [mean(M)]; markersize=8.0)
        errorbars!(ax, [2.0], [mean(M)], [std(M)]; whiskerwidth=5.0, linewidth=2.0)
        xlims!(ax, 0.0, 3.0)
    end

    # Plot summary stats for mean rates
    map(enumerate(zip(axs_σ, titles, results))) do (idx, (ax, title, result))
        Σ = map(zip(result, colors)) do (r, c)
            σ = std(r[t .> 0.1])
            scatter!(ax, [1.0], [σ]; color=c)
            σ
        end
        scatter!(ax, [2.0], [mean(Σ)]; markersize=8.0)
        errorbars!(ax, [2.0], [mean(Σ)], [std(Σ)]; whiskerwidth=5.0, linewidth=2.0)
        xlims!(ax, 0.0, 3.0)
    end

    axs[1, 2].title = "Inst. rate KDE"
    axs[1, 3].title = "μ"
    axs[1, 4].title = "σ"
    colsize!(fig.layout, 2, Relative(0.20))
    colsize!(fig.layout, 3, Relative(0.10))
    colsize!(fig.layout, 4, Relative(0.10))
    ylims!.(axs[:, 1:3], ylims[1]...)
    ylims!.(axs[:, 4], ylims[2]...)
    neaten_grid!(axs[:, 1:3])
    hidexdecorations!.(axs[1:2, 4]; ticklabels=true, ticks=false, grid=false)
    fig
end

function doplot_μ_histograms(stage="hsr", xlim=[0.0, 175.0], usecache=true)
    # Collect data
    fn = "\\home\\daniel\\cl_cache\\08-03-2023_temp.jld2"
    if isfile(fn) & usecache
        results = load(fn)["data"]
    else
        (noavg, internal), external = test(0.5, 100e3, 500, stage)
        results = [noavg, external, internal]
        save(fn, Dict("data" => results))
    end

    ## Plot
    # Create variables we need
    results = [noavg, external, internal]
    titles = ["No averaging", "External averaging", "Internal averaging"]

    # Create figures and axs
    fig = Figure()
    axs = [Axis(fig[i, 1]) for i in eachindex(results)]

    # Plot histogram of mean rates
    map(enumerate(zip(axs, titles, results))) do (idx, (ax, title, result))
        M = map(result) do r
            mean(r[t .> 0.1])
        end
        ax.title = title
        hist!(ax, M; bins=LinRange(xlim..., 50))
    end
    xlims!.(axs, xlim...)
    neaten_grid!(axs)
    axs[end].xlabel = "Mean rate (sp/s)"

    fig
end


# Do plots
fig = doplot_waveform_examples("hsr")

fig = doplot_waveform_examples("lsr", ([0.0, 5.0], [0.0, 2.0]))

doplot_μ_histograms()