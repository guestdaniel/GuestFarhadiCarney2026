export load_Warren1989a_Fig5b, plot_Warren1989a_Fig5b,
       load_Warren1989a_Fig6, plot_Warren1989a_Fig6,
       query_cas_threshold_Warren1989a_Fig6, query_cas_threshold_env_Warren1989a_Fig6,
       load_Warren1989a_Fig7, plot_Warren1989a_Fig7,
       load_Warren1989b_Fig2, plot_Warren1989b_Fig2,
       query_cas_ΔL_Warren1989b_Fig2, query_cas_ΔL_env_Warren1989b_Fig2,
       load_Warren1989b_Fig6, plot_Warren1989b_Fig6
"""
    load_Warren1989a_Fig5b()

Loads data scraped using WebPlotDigitizer from Figure 5b of:

Warren, E. H., & Liberman, M. C. (1989). Effects of contralateral sound on auditory-nerve 
responses. I. Contributions of cochlear efferents. Hearing Research, 37(2), 89–104. 
https://doi.org/10.1016/0378-5955(89)90032-4

The data is returned in the form of a data frame with three columns:
- `level`: Sound level of the contralateral suppressor (dB SPL)
- `rate`: Firing rate of the ipsilateral nerve fiber (sp/s)
- `condition`: "probe" for ipsi-only or "supp" for ipsi w/ contra suppressor
"""
function load_Warren1989a_Fig5b()
    # Loop over fiber types and load each data set separately
    dfs = map(["probe", "supp"]) do condition
        # Load file
        fn = projectdir("data", "Warren1989a_Fig5b_$(condition).csv")
        df = DataFrame(CSV.File(fn; header=false))
        rename!(df, :Column1 => :level, :Column2 => :rate)
        df[!, :condition] .= condition
        df
    end
    return vcat(dfs...)
end

function plot_Warren1989a_Fig5b(
    df=load_Warren1989a_Fig5b(); 
    size=(600, 400),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    plot_legend=false,
    linewidth=3.0,
)
    # Plot each dataset
    map(zip(["probe", "supp"], [:solid, :dash])) do (condition, linestyle)
        # Subset data
        df_sub = @subset(df, :condition .== condition)
        x = df_sub.level
        y = df_sub.rate

        # Plot
        lines!(ax, x, y; label=condition, linewidth=linewidth, linestyle=linestyle, color=:black)
    end

#    vlines!(ax, 40.0; color=:red)

    # Add legend
    if plot_legend
        axislegend(; position=:lt)
    end

    # Add ticks and labels
    ax.xticks = 0:20:100
    ax.xlabel = "Elicitor level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"

    # Adjust limits
    xlims!(ax, 0.0, 105.0)
    ylims!(ax, 0.0, 250.0)

    # Display and return
    display(fig)
    return fig, ax
end

"""
    load_Warren1989a_Fig6()

Loads data scraped using WebPlotDigitizer from Figure 6 of:

Warren, E. H., & Liberman, M. C. (1989). Effects of contralateral sound on auditory-nerve 
responses. I. Contributions of cochlear efferents. Hearing Research, 37(2), 89–104. 
https://doi.org/10.1016/0378-5955(89)90032-4

The data is in the form of a data frame with three columns:
- `cf`: Characteristic frequency of auditory-nerve fiber (kHz)
- `threshold`: Threshold for CAS suppression (dB SPL)
- `fiber_type`: "hsr", "msr", or "lsr"
"""
function load_Warren1989a_Fig6()
    # Loop over fiber types and load each data set separately
    dfs = map(["hsr", "msr", "lsr"]) do fiber_type
        # Load file
        fn = projectdir("data", "Warren_1989a_Fig6_$(uppercase(fiber_type)).csv")
        df = DataFrame(CSV.File(fn; header=false))
        rename!(df, :Column1 => :cf, :Column2 => :threshold)
        df[!, :fiber_type] .= fiber_type
        df
    end
    return vcat(dfs...)
end

function plot_Warren1989a_Fig6(
    df=load_Warren1989a_Fig6(); 
    size=(600, 400),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    normalized=false,
    markersize=14.0,
    linewidth=3.0,
    colors=colorschemes[:Dark2_8][[1, 2, 3]],
    plot_legend=false,
    plot_trendline=true,
    plot_envelope=true,
)
    # Create figure
    if !normalized
        xlims!(ax, 1.0, 2.0)
        ax.xscale=log10
        ax.xminorticksvisible=true
        ax.xminorticks=IntervalsBetween(9)
    end

    # Plot each dataset with scatter
    map(zip(["hsr", "msr", "lsr"], colors)) do (type, color)
        # Subset data
        df_sub = @subset(df, :fiber_type .== type)
        x = df_sub.cf
        y = df_sub.threshold

        # Plot
        scatter!(ax, x, y; label=type, color=(color, 0.5), markersize=markersize)
    end

    # Add trendline
    if plot_trendline
        x = df.cf
        y = df.threshold
        x̂, ŷ = smooth(log2.(x), y, quadratic)
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth)
    end
   
    if plot_envelope
        x = df.cf
        y = df.threshold
        x̂, ŷ = envelope(log2.(x), y, quadratic; ratio=200.0, offset=-6.0)  # envelope at 3 dB down
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth/2, linestyle=:dash)
    end

    # Add legend
    if plot_legend
        axislegend(; position=:lt)
    end

    # Add ticks and labels
    ax.xticks = ([0.1, 1.0, 10.0], ["0.1", "1", "10"])
    ax.yticks = 20.0:10.0:100.0
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "CAS threshold (dB SPL)"

    # Adjust limits
    xlims!(ax, 0.1, 30.0)
    ylims!(ax, 20.0, 90.0)

    # Display and return
    display(fig)
    return fig, ax
end

function query_cas_threshold_Warren1989a_Fig6(cf)
    df = load_Warren1989a_Fig6()
    x = df.cf
    y = df.threshold
    x̂, ŷ = smooth(log2.(x), y, quadratic)
    x̂ = 1e3(2 .^ x̂)
    ŷ[argmin(abs.(cf .- x̂))]
end

function query_cas_threshold_env_Warren1989a_Fig6(cf)
    df = load_Warren1989a_Fig6()
    x = df.cf
    y = df.threshold
    x̂, ŷ = envelope(log2.(x), y, quadratic)
    x̂ = 1e3(2 .^ x̂)
    ŷ[argmin(abs.(cf .- x̂))]
end

"""
    load_Warren1989a_Fig7()

Loads data scraped using WebPlotDigitizer from Figure 7 of:

Warren, E. H., & Liberman, M. C. (1989). Effects of contralateral sound on auditory-nerve 
responses. I. Contributions of cochlear efferents. Hearing Research, 37(2), 89–104. 
https://doi.org/10.1016/0378-5955(89)90032-4

The data is in the form of a data frame with three columns:
- `cf`: Characteristic frequency of auditory-nerve fiber (kHz)
- `magnitude`: Magnitude of CAS suppression (%)
- `fiber_type`: "hsr", "msr", or "lsr"
"""
function load_Warren1989a_Fig7()
    # Loop over fiber types and load each data set separately
    dfs = map(["hsr", "msr", "lsr"]) do fiber_type
        # Load file
        fn = projectdir("data", "Warren1989a_Fig7_$(uppercase(fiber_type)).csv")
        df = DataFrame(CSV.File(fn; header=false))
        rename!(df, :Column1 => :cf, :Column2 => :magnitude)
        df[!, :fiber_type] .= fiber_type
        df
    end
    return vcat(dfs...)
end

function plot_Warren1989a_Fig7(
    df=load_Warren1989a_Fig7(); 
    size=(500, 400),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    normalized=false,
    markersize=14.0,
    linewidth=3.0,
    colors=colorschemes[:Dark2_8][[1, 2, 3]],
    plot_legend=false,
    plot_trendline=true,
    plot_envelope=true,
)
    # Create figure
    if !normalized
        xlims!(ax, 1.0, 2.0)
        ax.xscale=log10
        ax.xminorticksvisible=true
        ax.xminorticks=IntervalsBetween(9)
    end

    # Add thick gray hline at zero
    hlines!(ax, [0.0]; color=:gray, linewidth=2.0)

    # Plot each dataset
    map(zip(["hsr", "msr", "lsr"], colors)) do (type, color)
        # Subset data
        df_sub = @subset(df, :fiber_type .== type)
        x = df_sub.cf
        y = df_sub.magnitude

        # Plot
        scatter!(ax, x, y; label=type, color=(color, 0.5), markersize=markersize)
    end

    # Add trendline
    if plot_trendline
        # Extract data
        x = df.cf
        y = df.magnitude
        x̂, ŷ = smooth(log2.(x), y, quadratic)
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth)
    end

    if plot_envelope
        # Extract data
        x = df.cf
        y = df.magnitude
        x̂, ŷ = envelope(log2.(x), y, quadratic; mode=:upper, offset=8.0, ratio=200.0)
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth/2, linestyle=:dash)
    end

    # Add legend
    if plot_legend
        axislegend(; position=:lt)
    end

    # Add ticks and labels
    ax.xticks = ([0.1, 1.0, 10.0], ["0.1", "1", "10"])
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "ΔR (% re: probe-alone)"

    # Adjust limits
    xlims!(ax, 0.1, 30.0)
    ylims!(ax, 100.0, -5.0) 

    # Display and return
    display(fig)
    return fig, ax
end

"""
    load_Warren1989b_Fig2()

Loads data scraped using WebPlotDigitizer from Figure 2 of:

Warren, E. H., and Liberman, M. C. (1989). “Effects of contralateral sound on auditory-nerve 
responses. II. Dependence on stimulus variables,” Hearing Research, 37, 105–121. 
doi:10.1016/0378-5955(89)90033-6

The data is returned in the form of a data frame with two columns:
- `cf`: CF of the afferent neuron (kHz)
- `ΔL`: Shift in RLF due to CAS (dB)
"""
function load_Warren1989b_Fig2()
    # Load file
    fn = projectdir("data", "Warren1989b_Fig2.csv")
    df = DataFrame(CSV.File(fn; header=false))
    rename!(df, :Column1 => :cf, :Column2 => :ΔL)
    df
end

function plot_Warren1989b_Fig2(
    df=load_Warren1989b_Fig2(); 
    size=(600, 400),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]),
    markersize=14.0,
    linewidth=3.0,
    plot_trendline=true,
    plot_envelope=true,
    alpha=0.5,
    color=:black,
)
    # Create figure
    xlims!(ax, 1.0, 2.0)
    ax.xscale=log10
    ax.xminorticksvisible=true
    ax.xminorticks=IntervalsBetween(9)

    # Plot data
    x = df.cf
    y = df.ΔL
    scatter!(ax, x, y; color=(color, alpha), markersize=markersize)

    # Add trendline
    if plot_trendline
        x̂, ŷ = smooth(log2.(x), y, quadratic; x₀=[log2(2.0), 1.0, 2.0])
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth)
    end

    if plot_envelope
        x̂, ŷ = envelope(log2.(x), y, quadratic; mode=:upper, ratio=100.0, offset=3.0)  # Envelope +3dB
        lines!(ax, 2 .^ x̂, ŷ; color=:gray, linewidth=linewidth/2, linestyle=:dash)
    end

    # Add ticks and labels
    ax.xticks = ([0.1, 1.0, 10.0], ["0.1", "1", "10"])
    ax.xlabel = "CF (kHz)"
    ax.ylabel = "CAS ΔL (dB)"

    # Adjust limits
    xlims!(ax, 0.1, 30.0)
    ylims!(ax, 0.0, 17.0)

    # Display and return
    display(fig)
    return fig, ax
end

function query_cas_ΔL_Warren1989b_Fig2(cf)
    df = load_Warren1989b_Fig2()
    x = df.cf
    y = df.ΔL
    x̂, ŷ = smooth(log2.(x), y, quadratic)
    x̂ = 1e3(2 .^ x̂)
    ŷ[argmin(abs.(cf .- x̂))]
end

function query_cas_ΔL_env_Warren1989b_Fig2(cf)
    df = load_Warren1989b_Fig2()
    x = df.cf
    y = df.ΔL
    x̂, ŷ = envelope(log2.(x), y, quadratic; mode=:upper)
    x̂ = 1e3(2 .^ x̂)
    ŷ[argmin(abs.(cf .- x̂))]
end


"""
    load_Warren1989b_Fig6()

Loads data scraped using WebPlotDigitizer from Figure 6 of:

Warren, E. H., and Liberman, M. C. (1989). “Effects of contralateral sound on auditory-nerve 
responses. II. Dependence on stimulus variables,” Hearing Research, 37, 105–121. 
doi:10.1016/0378-5955(89)90033-6

The data is returned in the form of a data frame with three columns:
- `freq`: Frequency of the contralateral suppressor (kHz)
- `rate`: Firing rate of the ipsilateral nerve fiber (sp/s)
- `condition`: "probe" for ipsi-only or "supp" for ipsi w/ contra suppressor
"""
function load_Warren1989b_Fig6()
    # Loop over fiber types and load each data set separately
    dfs = map(["probe", "supp"]) do condition
        # Load file
        fn = projectdir("data", "Warren1989b_Fig6_$(condition).csv")
        df = DataFrame(CSV.File(fn; header=false))
        rename!(df, :Column1 => :freq, :Column2 => :rate)
        df[!, :condition] .= condition
        df
    end
    return vcat(dfs...)
end

function plot_Warren1989b_Fig6(
    df=load_Warren1989b_Fig6(); 
    size=(500, 400),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]; xscale=log10, xminorticksvisible=true, xminorticks=IntervalsBetween(9)),
    plot_legend=false,
    linewidth=2.0,
    markersize=5.0,
    colors=[:black, :red],
)
    # Plot each dataset
    map(zip(["probe", "supp"], [:solid, :dash], colors)) do (condition, linestyle, c)
        # Subset data
        df_sub = @subset(df, :condition .== condition)
        x = df_sub.freq
        y = df_sub.rate

        # Plot
       scatter!(ax, x, y; label=condition, color=c, markersize=markersize)
       x̂, ŷ = smooth(log.(x), y, :loess; span=0.5)
       lines!(ax, exp.(x̂), ŷ; linewidth=linewidth, color=c)
    end

    # Add legend
    if plot_legend
        axislegend(; position=:lt)
    end

    # Add ticks and labels
    ax.xticks = [0.1, 1.0, 10.0]
    ax.xlabel = "Elicitor frequency (kHz)"
    ax.ylabel = "Firing rate (sp/s)"

    # Adjust limits
    xlims!(ax, 0.1, 10.0)
    ylims!(ax, 0.0, 250.0)

    # Display and return
    display(fig)
    return fig, ax
end
