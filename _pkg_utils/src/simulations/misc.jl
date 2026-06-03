export plot_nf

"""
    plot_nf(x::Vector, y::Vector{<:Vector})

Plot neurometric function between IV x and DV y.

Note we assume that x starts with a standard stimulus level of the IV. We assume that y
is a vector of which each element is a vector of the DV value across multiple repetitions of the stimulus at that level of the IV.
"""
function plot_nf(
    x::Vector{Float64}, 
    y::Vector{Vector{Float64}};
    ylims=[],
    fig=Figure(),
    ax=Axis(fig[1, 1]),
    color=:black,
    linewidth=1.0,
)
    # If nanmode == "include", leave values in y intact; 
    # Compute summary statistics
    μ_unmod = mean(y[1])
    σ_unmod = std(y[1])
    μ = map(mean, y[2:end])
    σ = map(std, y[2:end])

    # Plot
    errorbars!(ax, [minimum(x[2:end]) - 5.0], [μ_unmod], [σ_unmod]; color=color, whiskerwidth=5.0)
    scatter!(ax, [minimum(x[2:end]) - 5.0], [μ_unmod]; color=color)
    hlines!(ax, [μ_unmod]; linestyle=:dash, color=color, linewidth=linewidth)
    errorbars!(ax, x[2:end], μ, σ; color=color, whiskerwidth=5.0)
    scatter!(ax, x[2:end], μ; color=color)

    # Set limits
    if ~isempty(ylims)
        ylims!(ax, ylims...)
    end

    # Set labels
    ax.xlabel = "Mod. depth (dB)"
    ax.ylabel = "DV value"

    # Display
    display(fig)

    # Render
    fig, ax
end

