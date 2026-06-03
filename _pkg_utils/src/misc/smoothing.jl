export gwo, quickgwo, loper, quicksmooth_xy, quickquad, smooth, envelope, cubic, quadratic, logistic, slidestat, linear, gaussian, peaknorm_gaussian, smoothparam, gridparametric, explainparams, logistic3, invlogistic3, logistic_good

# Define a variety of useful parametric forms and default parameter values for them (via the
# `initval` function)
function cubic(x, θ)
    @. θ[1] + θ[2] * x + θ[3] * x^2 + θ[4] * x^3
end

function quadratic(x, θ)
    @. θ[1] + θ[2] * x + θ[3] * x^2
end

function linear(x, θ)
    @. θ[1] + θ[2] * x
end

function logistic(x, θ)
    @. θ[3] * (1 / (1 + exp(-(θ[1] + x * θ[2])))) + θ[4]
end

function logistic3(x, θ)
    @. θ[3] * (1 / (1 + exp(-(θ[1] + x * θ[2]))))
end

function logistic_good(x, θ)
    @. (θ[4]-θ[3]) * (1 / (1 + exp(-(-θ[1] * θ[2] + x * θ[2])))) + θ[3]
end

function invlogistic3(y, θ)
    @. (-θ[1] - log(θ[3]/y - 1))/θ[2]
end

function peaknorm_gaussian(x, μ=0.0, σ=1.0, w=1.0, o=0.0) 
    w * exp(-(x-μ)^2/(2σ^2)) + o
end

function gaussian(x::T, μ::T=0.0, σ::T=1.0, w::T=1.0, o::T=0.0) where {T <: AbstractFloat}
    w * 1/sqrt(2π * σ^2) * exp(-(x-μ)^2/(2σ^2)) + o
end

function gaussian(x::T, θ::Vector{T}) where {T <: AbstractFloat}
    gaussian.(x, θ[1], θ[2], θ[3], θ[4])
end

function initval(func::Function)
    if func == linear
        [0.0, 1.0]
    elseif func == quadratic
        [0.0, 1.0, 1.0]
    elseif func == cubic
        [0.0, 1.0, 1.0, 1.0]
    elseif func == logistic
        [0.0, 1.0, 1.0, 0.0]
    elseif func == logistic3
        [0.0, 0.0, 1.0]
    elseif func == logistic_good
        [0.0, 1.0, 1.0, 0.0]
    elseif func == gaussian
        [0.0, 1.0, 1.0, 0.0]
    end
end

function explainparams(func::Function)
    if func == linear
        ["Intercept", "Slope"]
    elseif func == quadratic
        ["Intercept", "x term", "x^2 term"]
    elseif func == cubic
        [0.0, 1.0, 1.0, 1.0]
    elseif func == logistic
        ["Intercept", "Slope", "Scale", "Offset"]
    elseif func == logistic3
        ["Intercept", "Slope", "Scale"]
    elseif func == gaussian
        [0.0, 1.0, 1.0, 0.0]
    end
end

mse(x, y) = mean(sum((x .- y) .^2))



# Define a `smooth` function and an `envelope` function that easily smooth data or extract
# an envelope, returning an interpolated x̂ and ŷ
"""
    smooth(x, y[, method=:loess; kwargs...])

Quick smooth data, returning x̂ and ŷ.

Smoothing options include:
- `:loess` for local regression, with kwargs `span`
- `:quadratic` for quadratic fit, with params offset and weights for powers of x
- `:logistic` for logistic fit, with params offset, slope, and scaling factor
- `:free` for free-form smoothing, with kwarg `ratio` 
"""
function smooth(x, y, method::Function=cubic; kwargs...)
    i = sortperm(x)
    x_sort = x[i]
    y_sort = y[i]
    quickparametric(x_sort, y_sort, method; kwargs...)[2:3]
end

function smoothparam(x, y, method::Function=cubic, approach="optim", args...; kwargs...)
    i = sortperm(x)
    x_sort = x[i]
    y_sort = y[i]
    if approach == "optim"
        quickparametric(x_sort, y_sort, method; kwargs...)
    else
        gridparametric(x_sort, y_sort, method, args...; kwargs...)
    end
end

function smooth(x, y, method::Symbol; kwargs...)
    i = sortperm(x)
    x_sort = x[i]
    y_sort = y[i]
    if method == :loess
        quickloess(x_sort, y_sort; kwargs...)
    elseif method == :free
        quickfree(x_sort, y_sort; kwargs...)
    end
end


"""
    envelope(x, y[, method=:loess; kwargs...])

Quick lower or upper envelope extraction, returning x̂ and ê (predicted envelope)
"""
function envelope(x, y, method::Function=quadratic; kwargs...)
    # Sort data in ascending x order
    i = sortperm(x)
    x_sort = x[i]
    y_sort = y[i]

    # Dispatch to method
    quickparametric_env(x_sort, y_sort, method; kwargs...)
end


# Implement the machinery for `smooth` and `envelope` via a few generic functions and
# a few specific more specific ones (for loess and sliding window estimation)
"""
    quickparametric(x, y[, func::Function=quadratic; x₀=[], n_pt=100])

Estimates parametric fit between `x` and `y` given functional form `func`
"""
function quickparametric(x, y, func::Function=quadratic; x₀=[], n_pt=1000)
    x̂ = LinRange(extrema(x)..., n_pt)
    f(θ) = mse(y, func(x, θ))
    θ̂ = Optim.minimizer(optimize(f, isempty(x₀) ? initval(func) : x₀; autodiff=ADTypes.AutoForwardDiff(), method=Newton()))
    return θ̂, x̂, func(x̂, θ̂)
end

function quickloess(x, y; n_pt=100, kwargs...)
    m = loess(x, y; kwargs...)
    x_hat = LinRange(extrema(x)..., n_pt)
    return x_hat, predict(m, x_hat)
end

function quickfree(x, y; ratio=0.2, interp_factor=1)
    # # Choose x̂ and identify which x̂ are in x
    # δ = (maximum(x) - minimum(x))/(interp_factor*(length(x)-1))
    # x̂ = minimum(x):δ:maximum(x)
    # if interp_factor == 1
    #     inx = fill(true, length(x̂))
    # else
    #     inx = mod.(1:length(x̂), interp_factor) .== 1
    # end

    # # Start ŷ at closest points by nearest neighbors
    # ŷ = fill(mean(y), length(x̂))
    # nearest = Interpolations.interpolate((x,), y, Gridded(Constant()))
    # ŷ = nearest[x̂] .- 5.0

    # # Define optimization function
    # # f(θ) = mean(diff(diff(θ)) .^ 2) + ratio*mse(y, θ[inx])
    # f(θ) = mean(diff(θ) .^ 2)# + ratio*mse(y, θ[inx])

    # # Minimzie and return fit
    # θ̂ = Optim.minimizer(optimize(f, ŷ; autodiff=:forward, iterations=50000))
    # return x̂, θ̂
    @warn "bad!"
end

"""
    gridparametric(x, y[, func::Function=quadratic; x₀=[], n_pt=100])

Estimates parametric fit between `x` and `y` given functional form `func` using grid search
"""
function gridparametric(x, y, func::Function, p...; n_pt=100)
    x̂ = LinRange(extrema(x)..., n_pt)
    f(θ) = mse(y, func(x, θ))
    P = collect(Iterators.product(p...))
    idx_min = argmin(map(f, P))
    return P[idx_min], x̂, func(x̂, P[idx_min]) 
end

"""
    quickparametric_env(x, y[, func::Function=quadratic; x₀=[], n_pt=100])

Estimates parametric envelope of `y` given functional form `func`
"""
function quickparametric_env(
    x, 
    y, 
    func::Function=quadratic; 
    x₀=[], 
    n_pt=100, 
    mode=:lower, 
    ratio=50.0, 
    exponent=2.0, 
    offset=0.0,
)
    # Sample x-axis over observed extrema with n_pt sample points
    x̂ = LinRange(extrema(x)..., n_pt)

    # Define loss function
    if mode == :lower
        f = θ -> mse(y, func(x, θ)) + ratio * sum(max.(0.0, func(x, θ) .- offset .- y) .^ exponent)
    elseif mode == :upper
        f = θ -> mse(y, func(x, θ)) + ratio * sum(max.(0.0, y .- (func(x, θ) .- offset)) .^ exponent)
    end

    # Optimize parameters given initial guess x₀ or default based on `initval`
    θ̂ = Optim.minimizer(optimize(f, isempty(x₀) ? initval(func) : x₀; autodiff=:forward))

    # Return
    return x̂, func(x̂, θ̂)
end


"""
    slidestat(x, y[, stat; len=0.3, jump=len/50, cutoff=30, fill=false, fillval=NaN])

Compute a statistic `stat` over a sliding window of width `len` and jump `jump`

`stat` is a function that takes in a vector of y values and returns a scalar. The default
`stat` is `median`. This function computes `stat` in a sliding window that moves by `jump`
proportion of the range of `x` and is length `len` proportion of the range of `x`. If 
`fill` if true, windows in which there are fewer than `cutoff` data points will be filled
with `fillval`.
"""
function slidestat(
    x, 
    y, 
    stat::Function=median; 
    len=0.3, 
    jump=len/10, 
    cutoff=30, 
    fill=false, 
    fillval=NaN,
    n_boot=1,
)
    # Define jump and len distance in absolute x units
    jump_abs = jump * (maximum(x) - minimum(x))
    len_abs = len * (maximum(x) - minimum(x))
    n_step = Int(ceil((maximum(x) - minimum(x)) / jump_abs))

    # Create empty vectors that we'll push into
    x̂ = minimum(x):jump_abs:maximum(x)
    ŷ = zeros(n_step, n_boot)

    # Loop over bootstraps
    for idx_boot in 1:n_boot
        # Resample data
        idx_resamp = rand(1:length(y), length(y))
        x_resamp = x[idx_resamp]
        y_resamp = y[idx_resamp]

        # Set position back to minimum of x
        pos = minimum(x)

        # Loop over sliding window query points
        for idx_step in 1:n_step
            # Calculate bounds of window
            lo = pos - len_abs/2
            hi = pos + len_abs/2

            # Next, identify data points in window and warn if we have too few points 
            idxs = (x_resamp .> lo) .& (x_resamp .< hi)
            n_pt = sum(idxs)
            if n_pt < cutoff
                @warn "Less than 30 data points in window!"
            end

            # If we have too few points and fill is true, fill with fillval, otherwise calculate stat in window and add to data
            if (n_pt < cutoff) & fill
                ŷ[idx_step, idx_boot] = fillval
            else
                ŷ[idx_step, idx_boot] = stat(y_resamp[idxs])
            end

            # Increment pos by jump distance
            pos += jump_abs
        end
    end

    # Return x̂, ŷ, and σ̂
    return x̂, dropdims(mean(ŷ; dims=2); dims=2), dropdims(std(ŷ; dims=2); dims=2)
end