export explain_fft_scaling, synccoef_cont, d′, nhpp_thinning, dot_raster, calc_psth, plot_psth, axis_grid
export getimg, displayimg, tilecat, padcat, autolog10ticks, truncate_below, truncate_above, truncate_outside, d′_pop_unbiased

# Function to explain how fft scaling works in FFTW
function explain_fft_scaling()
    x = randn(10)
    y = fft(x)
    println("Let's take the signal: $x")
    println("x has RMS power $(rms(x))")
    println("x has sum of power $(sum(x .^ 2))")
    println("")

    println("*Parseval's theorem*")
    println("Parseval's theorem, in the digital domain, indicates that:")
    println("∑_n=1^N a[n] a*[n] = 1/N ∑_k=1^N A[k] A*[k]")
    println("Confirming this... we can see that LHS ∑_n=1^N a[n] a*[n] is $(sum(x .^ 2))")
    println("Likewise, RHS 1/N ∑_k=1^N A[k] A*[k] is $(1 / length(y) * sum(y .* conj(y)))")
    println("Confirming this, we can see that LHS is $(sum(x .^ 2)) is and the RHS is $(1 / length(y) * sum(y .* conj(y)))")
    println("")

    println("*Scaling of FFT and IFFT*")
    println("Next, we simply note that in Julia real.(ifft(fft(x))) .== x: $(real.(ifft(fft(x))) .≈ x)")
    println("")

    println("*Calculate time-domain level from spectrum?*")
    println("Suppose given spectrum y, how do we calculate the sound level from the spectrum s.t. level(ifft(y)) matches our calculation?")
    println("First, we can observe how to match RMS(x) and level(y)...")
    println("RMS from time domain $(rms(x)) == RMS from spectral domain $(real(sqrt(1/length(y) * mean(y .* conj(y)))))")
    println("Therefore...")
    rms_spec = real(sqrt(1 / length(y) * mean(y .* conj(y))))
    println("Level from time domain $(20*log10(rms(x) / 20e-6)) dB SPL == RMS from spectral domain $(20*log10(rms_spec/20e-6)) dB SPL")
end

function synccoef_cont(r::Vector{Float64}, f::Float64, fs::Float64; cutoff=5.0)
    # Calculate vector strength
    t = timevec(r, fs)
    vs = abs(1 / sum(r) * sum(r .* exp.(1im .* 2π .* f .* t)))

    # Return NaN if mean(r) < 5 sp/s, otherwise return vs
    if mean(r) < cutoff
        return NaN
    else
        return vs
    end
end

# Calculate d′ for independent samples from two normal distributions
function d′(y1::Vector, y2::Vector; nanmode="omit", fillval=NaN, n_min=20)
    # Optionally handle nans in y1 and y2
    if nanmode == "omit"
        y1 = skipnan(y1)
        y2 = skipnan(y2)
    elseif nanmode == "fill"
        y1[isnan.(y1)] .= fillval
        y2[isnan.(y2)] .= fillval
    end

    # Require that both y1 and y2 have at least n_min non-NaN elements, otherwise return NaN
    if (sum(.!isnan.(y1)) < n_min) | (sum(.!isnan.(y2)) < n_min)
        return NaN
    else
        μ₁ = mean(y1)
        μ₂ = mean(y2)
        σ₁ = std(y1)
        σ₂ = std(y2)
        σ_pooled = sqrt(0.5 * (σ₁^2 + σ₂^2))
        return (μ₂ - μ₁) / σ_pooled
    end
end

function d′(y1::Vector{<:Vector}, y2::Vector{<:Vector}; kwargs...)
    map(zip(y1, y2)) do (_y1, _y2)
        d′(_y1, _y2; kwargs...)
    end
end

# Calculate d′ for numerous different sets of samples assuming first element is standard
function d′(y::Vector{<:Vector}; kwargs...)
    map(x -> d′(y[1], x; kwargs...), y[2:end])
end

"""
    d′_pop_unbiased(s::Vector{Vector{Float64}}, t::Vector{Vector{Float64}}; nanmode="omit", fillval=NaN, mask=[])

Estimate population-level d' from single-channel data using unbiased estimator of optimal 
combination rule

Assumes that s and t are vectors of vectors, where the outer vector is over channels and the 
inner vector is over reps
"""
function d′_pop_unbiased(
    s::Vector{Vector{Float64}}, t::Vector{Vector{Float64}};
    nanmode="omit",
    fillval=NaN,
    mask=[],
)
    # Handle what to do about NaNs (if nanmode=="include", continue, otherwise handle)
    # if nanmode == "omit"
    #     s = map(x -> map(_x -> _x[.!isnan.(_x)], x), s)
    #     t = map(x -> map(_x -> _x[.!isnan.(_x)], x), t)
    # elseif nanmode == "fill"
    #     error("fill not implemented")
    # end

    # Estimate d' variance from standard responses, separately for each channel
    var_by_chan = map(s) do chan_s
        # Loop over and do boostrapping to estimate d' variance
        # We know that E(X) = 0 here and we wish to estimate Var(X)
        sens_boot = map(1:1000) do _
            y1 = StatsBase.sample(chan_s, length(chan_s); replace=true)
            y2 = StatsBase.sample(chan_s, length(chan_s); replace=true)
            d′(y1, y2)
        end

        # Return estimate of variance
        var(sens_boot)
    end

    # Calculate d' using standard formula for every CF and depth (n_chan x 1 vector)
    sens = map(x -> d′(x[1], x[2]), zip(s, t))

    # Pool across zero-filled d' at each depth according to optimal combination rule
    # Each x^2 estimator is replaced with x^2 - Var(x)

    # First, we find the CFs that are valid for inclusion: must have estimated variance 
    # for the standard boostrapped d′ (no NaN) and must have non-nan d′ point estimate
    # at every depth 
    if isempty(mask)
        idxs_incl = (.!isnan.(var_by_chan)) .& .!isnan.(sens)
    else
        idxs_incl = (.!isnan.(var_by_chan)) .& .!isnan.(sens) .& mask
    end
    sens_ss = sens[idxs_incl]
    var_by_chan_ss = var_by_chan[idxs_incl]

    # Now, we implement out semi-adhoc rule (subtract variance from each squared d' estimate
    # and take the square root, unless the input is negative, in which case we set the 
    # result to zero)
    sqrt(max(0.0, sum(sens_ss .^ 2 .- var_by_chan_ss)))
end

"""
    nhpp_thinning(λ)

Generate draws from nonhomogeneous Poisson process via thinning

Algorithm shamelessly copied from https://stats.stackexchange.com/questions/369288/nonhomogeneous-poisson-process-simulation
"""
function nhpp_thinning(λ, fs)
    # Determine length of simulation
    T = length(λ) / fs

    # Determine the rate of the HPP from which we will thin as the maximum of the passed λ
    λ_max = maximum(λ)

    # Calculate the thinning probability, which is the ratio between λ and λ max
    p = λ ./ λ_max

    # Create an empty vector for spike times
    spiketimes = Float64[]

    # Initialize variables
    t = 0.0  # time of spike
    N = 0    # number of spikes

    # Loop and generate spikes
    while t < T
        # Generate first variable and update time
        U = rand(Uniform(0.0, 1.0))
        t = t - (log(U) / λ_max)

        if t >= (T - (1 / fs))
            # If we exceed T, return spiketimes
            return spiketimes
        else
            # If we don't exceed T, decide if we spike
            U₂ = rand(Uniform(0.0, 1.0))
            idx = sampleat(t, fs)
            if U₂ <= p[idx]
                N += 1
                push!(spiketimes, t)
            end
        end
    end
end

"""
    dot_raster(spks)

Very quick dot-raster plot with no frills
"""
function dot_raster(spks)
    fig = Figure()
    ax = Axis(fig[1, 1])
    for (idx, spk) in enumerate(spks)
        scatter!(ax, spk, fill(float(idx), length(spk)))
    end
    fig
end

"""
    calc_psth(spks, [T, binsize; mode])

From vector of vector of spike times re: stimulus onset, return PSTH
"""
function calc_psth(spks::Vector{Vector{Float64}}, T::Float64, binsize::Float64=1e-3; mode="instrate")
    # Count number of spike trains
    N = length(spks)

    # Pool all spike trains together
    spks = vcat(spks...)

    # Branch based on mode ("count", "instrate")
    if mode == "instrate"
        # If we want instantaneous rates, we normalize by binsize and by number of repeats
        psth = fit(StatsBase.Histogram, spks, 0.0:binsize:T)
        psth.weights .= Int64.(round.(psth.weights ./ binsize ./ N))
    elseif mode == "count"
        # If we want counts, we don't normalize by anything
        psth = fit(StatsBase.Histogram, spks, 0.0:binsize:T)
    end
    return psth
end

function calc_psth(spks::Vector{Vector{Vector{Float64}}}, args...; kwargs...)
    map(eachspk -> calc_psth(eachspk, args...; kwargs...), spks)
end


function plot_psth(
    psth::StatsBase.Histogram;
    strokewidth=0.0,
)
    # Create figure 
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Plot psth
    x = 0.5 .* (psth.edges[1][2:end] .+ psth.edges[1][1:(end-1)])
    barplot!(ax, x, psth.weights; gap=0.0, strokecolor=:black, color=:gray70, strokewidth=strokewidth)

    # Adjust labels
    ax.xlabel = "Time (s)"
    ax.ylabel = "PSTH bin count"

    # Return figure
    fig
end

plot_psth(x::Vector{<:Vector}, T=length(x) / 100e3, binsize=1e-3, args...; kwargs...) = plot_psth(calc_psth(x, args...); kwargs...)



"""
    axis_grid()

Given iterables and a figure object, return convenience variables for creating grid figure 
"""
function axis_grid(
    itr1,
    itr2;
    sz_h=150,
    sz_v=150,
    fig::Figure=Figure(; size=(sz_h * length(itr2), sz_v * length(itr1))),
    extra_xlabelsvisible=false,
    extra_ylablesvisible=false,
    extra_xticklabelsvisible=false,
    extra_yticklabelsvisible=false,
    formatter_row=x -> "x = $(string(x))",
    formatter_col=y -> "y = $(string(y))",
    labels_row=true,
    labels_col=true,
    kwargs...
)
    # Create empty arrays to store params and axs
    params = Array{Tuple{eltype(itr1),eltype(itr2)}}(undef, length(itr1), length(itr2))
    axs = Array{Axis}(undef, length(itr1), length(itr2))

    # Loop through and record parameters and create axis, passing kwargs to Axis
    for i in eachindex(itr1)
        for j in eachindex(itr2)
            params[i, j] = (itr1[i], itr2[j])
            axs[i, j] = Axis(
                fig[i, j];
                ylabelvisible=(!extra_ylablesvisible & (j > 1)) ? false : true,
                xlabelvisible=(!extra_xlabelsvisible & (i < length(itr1))) ? false : true,
                yticklabelsvisible=(!extra_yticklabelsvisible & (j > 1)) ? false : true,
                xticklabelsvisible=(!extra_xticklabelsvisible & (i < length(itr1))) ? false : true,
                kwargs...
            )
        end
    end

    # If labels, add labels at outer margins
    if labels_row
        for i in eachindex(itr1)
            Label(fig[i, length(itr2)+1], formatter_row(itr1[i]); tellheight=false)
        end
    end
    if labels_col
        for j in eachindex(itr2)
            Label(fig[0, j], formatter_col(itr2[j]); tellwidth=false)
        end
    end
    return params, fig, axs
end

function axis_grid(itr; orientation=:horizontal, formatter=x -> string(x), kwargs...)
    if orientation == :horizontal
        axis_grid([NaN], itr; formatter_col=formatter, labels_row=false, kwargs...)
    else
        axis_grid(itr, [NaN]; formatter_row=formatter, labels_col=false, kwargs...)
    end
end

function getimg(fig)
    Makie.update_state_before_display!(fig)
    permutedims(Makie.colorbuffer(Makie.get_scene(fig)))
end

function displayimg(
    img::Matrix{ARGB32},
    fig=Figure(; resolution=(size(img))),
    ax=Axis(fig[1, 1]; yreversed=true, aspect=DataAspect()),
)
    hidespines!(ax)
    hidedecorations!(ax)
    image!(ax, img)
    return fig
end

function tilecat(imgs)
    n = Int(ceil(sqrt(length(imgs))))
    chunks = map(Iterators.partition(imgs, n)) do chunk
        if length(chunk) < n
            vcat(chunk..., zeros(ARGB32, size(chunk[1])[1] * (n - length(chunk)), size(chunk[1])[2]))
        else
            vcat(chunk...)
        end
    end
    hcat(chunks...)
end

function tilecat(imgs::Matrix{Matrix{ARGB32}})
    rows = map(row -> hcat(row...), eachrow(imgs))
    vcat(rows...)
end

function padcat(img1, img2)
    img1 = copy(img1)
    img2 = copy(img2)
    s1 = size(img1, 2)
    s2 = size(img2, 2)
    if s1 < s2
        img1 = hcat(img1, zeros(ARGB32, size(img1, 1), s2 - s1))
    elseif s1 > s2
        img2 = hcat(img2, zeros(ARGB32, size(img2, 1), s1 - s2))
    end
    return vcat(img1, img2)
end

function displayimg(
    img::Matrix{Matrix{ARGB32}},
    fig=Figure(; resolution=size(img[1, 1])),
)
    n_row, n_col = size(img)
    for i in 1:n_row
        for j in 1:n_col
            displayimg(img[i, j], fig[i, j])
        end
    end
    fig
end

function autolog10lims(x; pad=1)
    # Find minimum value of x and position lowest tick at decade below minimum
    xmin = minimum(x)
    val_low = 10^(floor(log10(xmin)) - pad)

    # Same for highest above
    xmax = maximum(x)
    val_high = 10^(ceil(log10(xmax)) + pad)

    return (val_low, val_high)
end

function autolog10ticks(x; pad=1)
    # Get low and high values
    tick_low, tick_high = autolog10lims(x; pad=pad)
    return autolog10ticks(tick_low, tick_high)
end

function autolog10ticks(tick_low, tick_high)
    # Position ticks evenly on log scale from tick_low to tick_high
    ticks = 10 .^ (log10(tick_low):1:log10(tick_high))
    ticklabels = map(ticks) do tick
        n_digits = floor(log10(tick))
        n_digits = n_digits >= 0 ? 0 : Int(round(abs(n_digits)))
        fmt = FormatSpec("0.$(n_digits)g")
        pyfmt(fmt, tick)
    end
    ticks, ticklabels
end

function truncate_below(x, val)
    x[x.>val]
end

function truncate_above(x, val)
    x[x.<val]
end

function truncate_outside(x, val_low, val_high)
    x[(x.>val_low).&(x.<val_high)]
end

