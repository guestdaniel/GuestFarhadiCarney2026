export slidwin, neurogram

"""
    slidwin(x[; windowfunc=hamming, dur_window=0.01, fs=100e3])

Compute sliding-window average of `x` over time

Returns a vector of times (indicating the time at the first sample in the window) and 
a vector of corresponding average values of `x` at that time, computing by weighting
samples over `dur_window` with a windowfunc produced by `windowfunc(len_window)`, summating,
and dividing by the sum of the window alone.
"""
function slidwin(
    x::Vector{<:AbstractFloat};
    windowfunc=hamming,
    dur_window::Float64=4e-3,
    dur_skip::Float64=dur_window / 4,
    fs::Float64=100e3,
    mode="avg",
)
    # Convert skip and window durations into lengths in samples
    len_window = samples(dur_window, fs)
    len_skip = samples(dur_skip, fs)

    # Calculate number of resulting windows
    n_sample = length(x)
    n_window = Int(round((n_sample - (len_window - len_skip)) / len_skip) - 1)

    # Calculate indices for the first sample of each window
    idxs = 1:len_skip:(n_sample-(len_window-len_skip))

    # Preallocate
    output = zeros(n_window)
    w = windowfunc(len_window)
    norm_w = sum(w)

    # Loop through time and calculate
    for idx_step in 1:n_window
        incl = idxs[idx_step]:(idxs[idx_step]+len_window-1)
        if mode == "avg"
            output[idx_step] = sum(x[incl] .* w) / norm_w
        elseif mode == "rms"
            output[idx_step] = rms(x[incl] .* w) / norm_w
        else
            output[idx_step] = mode(x[incl])
        end
    end

    # Return
    return timeat.(idxs[1:n_window], fs), output
end

function slidwin(x::Matrix{T}; kwargs...) where {T<:AbstractFloat}
    out = mapslices(z -> slidwin(z; kwargs...)[2], x; dims=1)
    slidwin(x[:, 1]; kwargs...)[1], out
end

function slidwin(x::Vector{Vector{T}}; kwargs...) where {T<:AbstractFloat}
    out = map(z -> slidwin(z; kwargs...), x)
    out[1][1], getindex.(out, 2)
end

"""
    neurogram(x)

Plot contents of x as colorplot
"""
function neurogram!(
    ax::Axis,
    cf::Vector{Float64},
    x::Matrix{Float64};
    windowfunc=hamming,
    dur_window::Float64=4e-3,
    fs::Float64=100e3,
    colorrange=(0.0, Inf),
    ratio_autoy=1.1,
    colormap=colorschemes[:viridis],
    colorscale=identity,
    yticks=[0.5, 1.0, 2.0, 4.0, 8.0, 16.0],
    kwargs...
)
    # Reduce from full data to smoothed/subsampled data
    t_smooth, x_smooth = slidwin(x; windowfunc=windowfunc, dur_window=dur_window, fs=fs)

    # Plot figure
    if colorrange == (0.0, Inf)
        colorrange = (0.0, ratio_autoy * maximum(x_smooth))
    end
    hm = heatmap!(ax, t_smooth, cf ./ 1e3, x_smooth; colorrange=colorrange, colormap=colormap, colorscale=colorscale)

    # Adjust labels and ticks
    ax.xlabel = "Time (s)"
    ax.ylabel = "CF (kHz)"
    ax.yticks = yticks
    hm
end

neurogram!(ax, cf, x::Vector{Vector{Float64}}; kwargs...) = neurogram!(ax, cf, hcat(x...); kwargs...)

function neurogram(args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1]; yscale=log10)
    hm = neurogram!(ax, args...; kwargs...)
    Colorbar(fig[1, 2], hm; label="Firing rate (sp/s)")
    fig
end