export saveplot,
       getie,
       fetch_model,
       fetch_model_baseline,
       plot_ie,
       MOCParams3,
       MOCParams4,
       MOCParams5,
       fetch_params,
       fetch_attr,
       hsr_threshold,
       lsr_threshold,
       hsr_midpoint,
       lsr_midpoint,
       cs_level,
       cohc_to_ΔL
export linescatter!, list, @parallel, mse, middle, vs, colormap_freq_bipolar,
    quicklook_gain_curves_pt, tradlogxax, quicksmooth, moc_nonlinearity, tradlogdoubleax

# Convenience function to save figures at 600 dpi
saveplot(name, fig) = save(projectdir("figs", name), fig; px_per_unit=6) 


# Function to generate "innervation exponent" at each CF location based on CF in Hz
# _, _, θ = fetch_fit_Liberman1990_Fig14a([log2(0.05), log2(100)])  # fit between CF correlate and OHC terminal area
θ = [2.003, 2.282, 296.820, -3.557]
getie(cf; shift=-1.5) = gaussian.(log2(cf/1e3 * 2^-shift), θ...) / gaussian(θ[1], θ...)


# Structure to encode key MOC model params, as easy alterantive format to provide a method
# to the GFC2024_FakeBin constructors
abstract type MOCParams end

@with_kw struct MOCParams3 <: MOCParams
    moc_weight_func::Function=x -> 5.0  # function that maps from CF (Hz) to weight
    moc_width::Float64=1.0
    moc_offset::Float64=0.0
    moc_beta::Float64=0.015
end

@with_kw struct MOCParams4 <: MOCParams
    moc_weight::Function=x -> 1.0
    moc_offset::Function=x -> 5.0
    moc_beta::Function=x -> 0.2
    moc_width::Float64=1.0
end

function Utilities.viz(params::MOCParams3; fig=Figure(), ax=Axis(fig[1, 1]), cf=2e3)
    # Plot rate vs moc nonlinearity output
    x = LinRange(0.0, 50.0, 1000)
    y = moc_nonlinearity.(x; w=params.moc_weight_func(cf), β=params.moc_beta, θ=params.moc_offset, minval=0.0, maxval=1.0)
    lines!(ax, x, y)

    # Add labels and such
    ax.xlabel = "Firing rate (sp/s)"
    ax.ylabel = "Gain factor"

    # Limits
    xlims!(ax, 0.0, maximum(x))
    ylims!(ax, 0.0, 1.05)

    # Return
    display(fig)
    fig, ax
end

function Utilities.viz(params::MOCParams4; fig=Figure(; size=(600, 200)), axs=[Axis(fig[1, 1]), tradlogxax(fig[1, 2]), tradlogxax(fig[1, 3])], cf=2e3)
    # Plot rate vs moc nonlinearity output
    x = LinRange(0.0, 50.0, 1000)
    y = moc_nonlinearity.(x; w=params.moc_weight(cf), β=params.moc_beta(cf), θ=params.moc_offset(cf), minval=0.0, maxval=1.0)
    lines!(axs[1], x, y)
    axs[1].xlabel = "Firing rate (sp/s)"
    axs[1].ylabel = "Gain factor"
    xlims!(axs[1], 0.0, maximum(x))
    ylims!(axs[1], 0.0, 1.05)

    # Plot beta vs CF
    x = LogRange(0.25, 16e3, 500)
    lines!(axs[2], x, params.moc_beta.(x .* 1e3))
    axs[2].xlabel = "CF (kHz)"
    axs[2].ylabel = "β"
    xlims!(axs[2], 0.1, 30.0)

    # Plot offset vs CF
    x = LogRange(0.25, 16e3, 500)
    lines!(axs[3], x, params.moc_offset.(x .* 1e3))
    axs[3].xlabel = "CF (kHz)"
    axs[3].ylabel = "Offset"
    xlims!(axs[3], 0.1, 30.0)

    # Return
    display(fig)
    fig, axs[1]
end

function Utilities.viz(params::Vector{MOCParams3})
    fig = Figure()
    ax = Axis(fig[1, 1])
    for param in params
        viz(param; fig=fig, ax=ax)
    end
    display(fig)
    return fig, ax
end

function Utilities.viz(params::Vector{MOCParams4})
    fig = Figure()
    ax = Axis(fig[1, 1])
    for param in params
        viz(param; fig=fig, ax=ax)
    end
    display(fig)
    return fig, ax
end

"""
    fetch_params(id::String="default")

Returns set of MOCParams3 based on ID string
"""
function fetch_params(id::String="default")
    @match id begin
        # "default" matches parameter defaults above
        "default" => MOCParams4()
        # "current" is current best overall result
        "current" => MOCParams4(;
            moc_weight=cf -> cf < 2e3 ? 5.0 : 5.0 * peaknorm_gaussian(log2(cf/2e3), 0.0, 2.0),
            moc_offset=cf -> 0.0,
            moc_beta=cf -> 0.015,
            moc_width=1.0,
        )
        # Best tries to compormise between threshold and magnitude
        "best" => MOCParams4(;
            moc_weight=cf -> cf < 2e3 ? 40.0 : 40.0 * peaknorm_gaussian(log2(cf/2e3), 0.0, 2.0),
            moc_offset=cf -> 0.0,
            moc_beta=cf -> 0.0016,
            moc_width=1.0,
        )
    end
end


# Functions to generate models using standard parameter choices
"""
    fetch_model(cf::Float64=1e3 [; stage="hsr", guardrail_mode="none", moc_weight_wdr=6.0])
    fetch_model(params::MOCParam3, cf::Float64=1e3 [; stage="hsr", guardrail_mode="none", moc_weight_wdr=6.0])

Generates a GFC2024_FakeBin object that uses standard parameters for the 2025 manuscript.
"""
function fetch_model(cf::Float64=1e3; kwargs...)
    fetch_model(
        fetch_params("current"),
        cf; 
        kwargs...
    )
end

function fetch_model(
    params::MOCParams3, 
    cf::Float64=1e3;
    stage="hsr",
    species="cat",
    n_cf=41,
    cf_range=3.0,  # octaves from lowest to highest CF around center CF
    cfs=round.(LogRange(cf*2^-(cf_range/2), cf*2^(cf_range/2), n_cf)),
    kwargs...
)
    GFC2024_FakeBin(
        cfs; 
        cf_ipsi=cf, 
        stage=stage, 
        moc_weight=params.moc_weight_func.(cfs), 
        moc_width=params.moc_width, 
        moc_offset=params.moc_offset, 
        moc_beta=params.moc_beta,
        species=species,
        kwargs...
    )
end

function fetch_model(
    params::MOCParams4, 
    cf::Float64=1e3;
    stage="hsr",
    species="cat",
    n_cf=41,
    cf_range=3.0,  # octaves from lowest to highest CF around center CF
    cfs=round.(LogRange(cf*2^-(cf_range/2), cf*2^(cf_range/2), n_cf)),
    kwargs...
)
    GFC2024_FakeBin(
        cfs; 
        cf_ipsi=cf, 
        stage=stage, 
        moc_weight=params.moc_weight.(cfs), 
        moc_width=params.moc_width, 
        moc_offset=params.moc_offset.(cfs), 
        moc_beta=params.moc_beta.(cfs),
        species=species,
        kwargs...
    )
end


"""
    fetch_model_baseline(cf; kwargs...)

Generates a GFC2024 object that uses standard parameters for the 2025 manuscript.
"""
function fetch_model_baseline(
    cf::Float64=1e3;
    stage="hsr",
    species="cat",
    cohc=1.0,
    kwargs...
)
    # Branch based on request for model type; 
    GFC2024(; 
        cf=[cf], 
        stage=stage, 
        moc_weight=[0.0],
        moc_width=0.0,
        cohc=[cohc],
        species=species,
        kwargs...
    )
end

function fetch_model_baseline(
    cf::Vector{Float64};
    stage="hsr",
    species="cat",
    cohc=[1.0],
    mode="multichannel",
    kwargs...
)
    GFC2024(; 
        cf=cf, 
        stage=stage, 
        moc_weight=zeros(size(cf)),
        moc_width=0.0,
        cohc=cohc,
        species=species,
        coi=mode=="singlechannel" ? middle(1:length(cf)) : 1:length(cf),
        kwargs...
    )
end

function fetch_model_baseline(model::Model; kwargs...)
    fetch_model_baseline(
        getcf(model); 
        stage=ipsi(model).stage, 
        species=ipsi(model).species, 
        cohc=ipsi(model).cohc[1],
        kwargs...
    )
end

# Functions to quickly query HSR and LSR midpoints and thresholds using reasonable
# level resolution (1 dB) and suitable level ranges. For visual support see:
#   genfig_supp_RLFs_debug(; config=config, fiber_type="hsr")
#   genfig_supp_RLFs_debug(; config=config, fiber_type="lsr")
function hsr_threshold(cf; config=Config())
    model = fetch_model_baseline(cf; stage="hsr")
    sim = RLF(model; levels=-5.0:1.0:45.0)
    Utilities.threshold(sim; config=config)
end

function lsr_threshold(cf; config=Config())
    model = fetch_model_baseline(cf; stage="lsr")
    sim = RLF(model; levels=-5.0:1.0:45.0)
    Utilities.threshold(sim; config=config)
end

function hsr_midpoint(cf; config=Config())
    model = fetch_model_baseline(cf; stage="hsr")
    sim = RLF(model; levels=-5.0:1.0:45.0)
    Utilities.midpoint(sim; config=config)
end

function lsr_midpoint(cf; config=Config())
    model = fetch_model_baseline(cf; stage="lsr")
    sim = RLF(model; levels=-5.0:1.0:45.0)
    Utilities.midpoint(sim; config=config)
end


"""
    fetch_attr(attr::Symbol, id::Symbol)

Fetches attribute (e.g., :linestyle) based on id (e.g., :model)
"""
function fetch_attr(attr::Symbol, id::Symbol)
    @match attr begin
        :color => @match id begin
            :model => :red
            :physio => :gray
            :lsr => ColorSchemes.Dark2_8[3]
            _ => :black
        end
    end
end

cs_level(level; level_min=0.0, level_max=100.0) = level == -Inf ? :black : get(colorschemes[:glasgow], (level-level_min)/level_max)

"""
    cohc_to_ΔL(vq; config=Config(), cf=1e3, species="human")

Transforms a query COHC value into an empirical ΔL 
"""
function cohc_to_ΔL(vq; config=Config(), cf=1e3, species="human")
    cohcs = LogRange(0.01, 1.0, 21)
    rlfs = map(cohcs) do cohc
        m = GFC2024(
            [cf]; 
            stage="hsr", 
            moc_weight=[0.0],
            moc_offset=[1000.0],
            species=species,
            cohc=[cohc],
        )
        RLF(m; levels=0.0:2.0:80.0)
    end
    ΔL = map(rlfs[1:(end-1)]) do rlf
        midpoint_shift(rlfs[end], rlf; config=config)
    end
    itp = linear_interpolation(cohcs, vcat(ΔL..., 0.0), extrapolation_bc=Line())
    itp(vq)
end

function linescatter!(
    ax,
    x,
    y;
    linewidth=1.0,
    markersize=14.0,
    linestyle=:solid,
    filled=true,
    label="",
    kwargs...
)
    lines!(ax, x, y; linewidth=linewidth, linestyle=linestyle, kwargs...)
    if filled
        scatter!(ax, x, y; markersize=markersize, label=label, kwargs...)
    else
        scatter!(ax, x, y; markersize=markersize, label=label, kwargs...)
        scatter!(ax, x, y; markersize=markersize / 3, kwargs..., color=:white)
    end
end

macro parallel(n=4)
    quote
        using Distributed
        if nprocs() == nworkers() == 1
            addprocs($n)
        elseif nworkers() < $n
            addprocs($n - nworkers())
        end
        @everywhere using Pkg
        @everywhere Pkg.activate(".")
        @everywhere using GuestFarhadiCarney2026
        @everywhere using Utilities
        @everywhere using Interpolations
        @info "Parallel pool established with $(nworkers()) workers, GuestFarhadiCarney2026 loaded!"
    end
end

mse(y, ŷ) = mean((y .- ŷ) .^ 2)

middle(x) = x[Int(ceil(length(x) / 2))]

vs(x, f=100.0, fs=100e3) = abs(1 / sum(x) * sum(x .* exp.(1im .* 2π .* f .* timevec(x, fs))))

function colormap_freq_bipolar(freqs, freq_ref; extent=1.0, gain=1.5, colormap=:berlin)
    # Map freqs into signed octave distance from freq_ref, capped by extent (i.e., if distance
    # is -1.5 octaves and extent is 1.0, distance should be capped to -1.0).
    dist = log2.(freqs ./ freq_ref)
    dist[abs.(dist).>extent] .= sign.(dist[abs.(dist).>extent]) .* extent

    # Map to colormap
    get(colorschemes[colormap], gain .* (dist .+ 1) ./ 2)
end

function tradlogxax(args...; kwargs...)
    Axis(
        args...;
        xscale=log10,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        xticks=[0.1, 1.0, 10.0],
        kwargs...
    )
end

function tradlogdoubleax(args...; kwargs...)
    Axis(
        args...;
        xscale=log10,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        xticks=[0.1, 1.0, 10.0],
        yscale=log10,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
        yticks=[0.1, 1.0, 10.0],
        kwargs...
    )
end


function quicklook_gain_curves_pt(
    config::Config=Config();
    freq=2e3,
    level=65.0,
    widths=[1.0, 2.0, 4.0],
    weights=[5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0],
)
    # Create pure-tone stimulus
    stim = PureTone(; freq=freq, dur=1.0, dur_ramp=2.5e-3, level=level)
    cf = LogRange(0.2e3, 20e3, 41)

    # Compute curves, mapping over product of weights and widths
    curves = pmap(Iterators.product(weights, widths)) do (weight, width)
        model = GFC2024(cf; moc_weight_ic=0.0, moc_width_wdr=width, moc_weight_wdr=weight, guardrail_mode="none")
        resp = Utilities._compute(model, stim)
        map(x -> mean(x[50000:100000]), resp["gainpostmix"])
    end

    # Plot curves
    css = [
        colorschemes[:Reds],
        colorschemes[:Blues],
        colorschemes[:Greens],
        colorschemes[:Purples],
        colorschemes[:Oranges],
    ]
    fig = Figure(; size=(350, 800))
    map(enumerate(zip(eachcol(curves), css))) do (idx, (curves_by_weight, cs))
        ax = Axis(fig[idx, 1]; xscale=log10, yscale=log10)
        colors = get(cs, LinRange(0.0, 1.0, length(curves_by_weight)))
        map(zip(curves_by_weight, colors)) do (curve, color)
            lines!(cf, curve; color=color)
        end
        ax.xticks = [0.5e3, 1e3, 2e3, 4e3, 8e3]
        ylims!(ax, 5e-2, 1.0)
        ax.yticks = [0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
        if idx < length(widths)
            hidexdecorations!(ax, ticks=false, grid=false)
        end
    end
    fig
end


"""
    moc_nonlinearity(x[; w=1.0, β=0.015, θ=40.0, minval=0.1, maxval=1.0])
    moc_nonlinearity(x, params, cf)

Compute output of pointwise MOC input-output nonlinearity

Computes output of rational input-output nonlinearity between LSR rates and output
gain factor, given weight on input `w`, slope factor `β`, threshold w.r.t. to weighted rate 
`θ`.
"""
function moc_nonlinearity(x; w=1.0, β=0.015, θ=40.0, minval=0.1, maxval=1.0)
    if (w * x) < θ
        return maxval
    else
        return ((maxval - minval) / (1.0 + (β * ((w * x) - θ))^2)) + minval
    end
end

function moc_nonlinearity(x, params, cf)
    # Extract parameters
    w = params.moc_weight(cf)
    β = params.moc_beta(cf)
    offset = params.moc_offset(cf)

    # Return
    moc_nonlinearity(x; w=w, β=β, θ=offset, minval=0.0, maxval=1.0)
end

