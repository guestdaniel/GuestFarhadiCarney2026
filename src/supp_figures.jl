export genfig_supp_cf_density,
       genfig_supp_cas_thresholds_parameters,
       genfig_supp_cas_thresholds_parameters_pcolor,
       genfig_supp_cas_single_curve_grid,
       genfig_supp_cas_rlf_shift_grid,
       genfig_supp_cas_ΔL_vs_level,
       genfig_supp_cas_ΔL_vs_level_condensed,
       genfig_supp_cas_ΔL_vs_level_β,
       genfig_survey_moc_nonlin,
       genfig_supp_elicitor_vs_probe_response_dynamics,
       genfig_supp_cas_suppression_growth,
       genfig_supp_level_growth_panel,
       genfig_supp_rlf_steepness,
       genfig_supp_ΔL_vs_threshold_params,
       genfig_supp_ΔL_vs_threshold_grid,
       genfig_supp_tc_diagnose,
       genfig_supp_c1_compression,
       genfig_supp_RLFs_debug


"""
    genfig_supp_cf_density()

Supplemental figure to analyze the effect of CF density on simulation results
"""
function genfig_supp_cf_density(; config=Config())
    # Params
    nos = [5, 11, 21, 41]
    colors = colorschemes[:Dark2_8]

    # Rate-level functions
    sims = map(nos) do no
        model = fetch_model(; n_cf=no)
        RLF(model, ProbePureToneElicitor; freq_elicitor=getcf(model)*1.1) 
    end
    r = map(sims) do sim
        @memo config simulate(sim)
    end
    fig = Figure()
    ax = Axis(fig[1, 1])
    map(zip(r, colors)) do (_r, c)
        scatter!(ax, level(sims[1]), _r; color=c, markersize=4.0)
        lines!(ax, level(sims[1]), _r; color=c)
    end
    ylims!(ax, 0.0, 300.0)
    ax.xlabel = "Elicitor level (dB SPL)"
    ax.ylabel = "Firing rate (sp/s)"
    fig
end


"""
    genfig_supp_cas_thresholds_parameters()

Plot CAS thresholds as function of model parameters (weight and IO offset [threshold])

Demonstrates, at a single example CF, that both the weight parameter and the offset
parameter contribute to setting CAS threshold. See also the
`genfig_supp_cas_thresholds_parameters_pcolor` for a colorplot showing the difference
between the nominal "target" CAS threshold (based on physio data) and the simulated CAS
thresholds. The colorplot makes it clear that there is a contiguous area for the two
parameters that achieve accurate CAS thresholds; then disambiguating these parameters may be
by achieved by comparing suprathreshold behavior. Comparing these plots at low and high CFs
also makes some other issues manifest, such as the relatively low CAS thresholds in the
model at high CFs compared to physiological data.
"""
function genfig_supp_cas_thresholds_parameters(; 
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = 2.0:1.0:20.0
    moc_offsets_factor = 0.0:1.0:8.0

    # Do sims
    curves = map(moc_offsets_factor) do offset
        map(moc_weights) do weight
            model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
            level_probe = midpoint(RLF(ipsi(model)); config=config)
            sim = SLF(model; levels_elicitor=25.0:1.0:60.0, level_probe=level_probe)
            threshold(sim; config=config)
        end
    end

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Add line indicating physiological threshold
    th = query_cas_threshold_Warren1989a_Fig6(cf)
    hlines!(ax, [th]; color=:gray)

    # Loop over different curves and plot
    lns = map(zip(curves, get(ColorSchemes.Purples, LinRange(0.2, 1.0, length(curves))))) do (curve, c)
        # Markers that are bad (NaN) should be handled separately
        idxs_nan = isnan.(curve)
        scatter!(ax, moc_weights[idxs_nan], 65.0 .* ones(sum(idxs_nan)); color=c, marker='↑', markersize=10.0)

        # Plot line 
        scatter!(ax, moc_weights, curve; color=c, markersize=7.0)
        lines!(
            ax, 
            quicksmooth_xy(moc_weights[.!(idxs_nan)], curve[.!(idxs_nan)]; span=0.4)...; 
            color=c
        )
    end
    Legend(fig[1, 2], lns, string.(moc_offsets_factor), "offset factor")
    
    # Add labels
    ax.xlabel = "MOC weight (param a.u.)"
    ax.ylabel = "CAS threshold (dB SPL)"

    # Adjust
    ylims!(ax, 25.0, 70.0)
    xlims!(ax, 0.0, 20.0)

    # Display and save
    display(fig)
end


function genfig_supp_cas_thresholds_parameters_pcolor(; 
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = 2.0:1.0:20.0
    moc_offsets_factor = 0.0:1.0:8.0

    # Do sims
    Θ = map(Iterators.product(moc_offsets_factor, moc_weights)) do (offset, weight)
        model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
        level_probe = midpoint(RLF(ipsi(model)); config=config)
        sim = SLF(model; levels_elicitor=25.0:1.0:60.0, level_probe=level_probe)
        threshold(sim; config=config)
    end

    # Determine target threshold at this CF; colormap is adjusted in comparison
    target = query_cas_threshold_Warren1989a_Fig6(cf)
    Θ = Θ .- target

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1])
    hm = heatmap!(
        ax, 
        moc_offsets_factor, 
        moc_weights, 
        Θ; 
        colorrange=(-30.0, 30.0), 
        colormap=ColorSchemes.RdBu,
    )
    Colorbar(fig[1, 2], hm; label="Sim. thr − physio. thr (dB)")

    # Add labels
    ax.xlabel = "IO nonlinearity offset (× weight)"
    ax.ylabel = "MOC weight"

    # Display and save
    display(fig)
end


"""
    genfig_supp_cas_single_curve_grid()

Plot CAS thresholds as function of model parameters.

Figure showing rudimentary visual "grid search" for appropriate efferent model parameters
based on resulting changes in CAS thresholds across CF range.
"""
function genfig_supp_cas_single_curve_grid(;
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = [5.0, 10.0, 20.0]
    moc_offsets_factor = [0.0, 2.0, 4.0, 8.0]
    fig = Figure(; size=(200*length(moc_offsets_factor), 200*length(moc_weights)))
    params, axs = axis_grid(
        moc_weights,          # rows
        moc_offsets_factor,   # cols
        formatter_row=x->"weight = $x",
        formatter_col=x->"offset factor = $x",
        fig
    )

    # Plot different options
    map(zip(params, axs)) do ((weight, offset), ax)
        model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
        sim = SLF(model)
        viz(sim; config=config, fig=fig, ax=ax)
        vlines!(ax, [threshold(sim; config=config)]; color=:black)
    end

    # Display and save
    display(fig)
end


"""
    genfig_supp_cas_rlf_shift_grid()

Plot shifted RLFs as function of model parameters, acoustic parameters, etc.
"""
function genfig_supp_cas_rlf_shift_grid(;
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = [5.0, 10.0, 20.0]
    moc_offsets_factor = [0.0, 2.0, 4.0, 8.0]
    fig = Figure(; size=(200*length(moc_offsets_factor), 200*length(moc_weights)))
    params, axs = axis_grid(
        moc_weights,          # rows
        moc_offsets_factor,   # cols
        formatter_row=x->"weight = $x",
        formatter_col=x->"offset factor = $x",
        fig
    )

    # Loop over parameters and generate
    map(zip(params, axs)) do ((weight, offset), ax)
        model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
        sim = RLF(model, ProbePureToneElicitor; level_elicitor=-Inf, level_spacing=5.0)
        sim2 = RLF(model, ProbePureToneElicitor; levels=level(sim))
        viz(sim; config=config, fig=fig, ax=ax)
        viz(sim2; config=config, fig=fig, ax=ax, color=:red)
    end

    # Display and save
    display(fig)
end

"""
    genfig_supp_cas_ΔL_vs_level()

Plot shifted RLFs as function of model parameters, acoustic parameters, etc.

This figure shows RLF shifts as a function of elicitor level for an on-CF elicitor, with 
different rows and columns indicating different combinations of weights and offsets. 
Generally, we see somewhat accelerating ΔL-elicitor level functions after CAS threshold is 
achieved; it is not completely clear if this shallow out at higher sound levels? To a first
approximation, the offset factor seems to mostly shift the curves rightward whereas the 
weight seems to increase curve slopes.
"""
function genfig_supp_cas_ΔL_vs_level(;
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = [5.0, 7.5, 10.0, 12.5]
    moc_offsets_factor = [0.0, 2.0, 4.0, 8.0]
    levels_elicitor = 20.0:5.0:80.0

    # Create figure and axes
    fig = Figure(; size=(200*length(moc_offsets_factor), 200*length(moc_weights)))
    params, axs = axis_grid(
        moc_weights,          # rows
        moc_offsets_factor,   # cols
        fig,                  # fig obj
        formatter_row=x->"weight = $x",
        formatter_col=x->"offset factor = $x",
        xlabel="Elicitor level (dB SPL)",
        ylabel="ΔL (dB)",
    )

    # Loop over parameters and generate
    curves = map(params) do (weight, offset)
        map(levels_elicitor) do level_elicitor
            model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
            sim1 = RLF(model, ProbePureToneElicitor; level_elicitor=-Inf, level_spacing=2.5)
            sim2 = RLF(model, ProbePureToneElicitor; level_elicitor=level_elicitor, levels=level(sim1))
            threshold_shift(sim1, sim2; config=config)
        end
    end

    # Plot
    map(zip(axs, curves)) do (ax, curve)
        scatter!(ax, levels_elicitor, curve; color=:black)
        lines!(ax, quickquad(levels_elicitor, curve)...; color=:black)
        ablines!(ax, [-20.0], [1.0]; color=:lightgray, linestyle=:dash)
        ablines!(ax, [-30.0, -40.0, -50.0, -60.0], [1.0, 1.0, 1.0, 1.0]; color=:lightgray, linestyle=:dot, linewidth=0.5)
    end

    # Adjust limits
    ylims!.(axs, 0.0, 35.0)

    # Display and save
    display(fig)
end

function genfig_supp_cas_ΔL_vs_level_condensed(;
    cf=2e3,
    config=Config(),
)  
    # Set up axes for weight sweep
    moc_weights = [5.0, 10.0, 15.0, 20.0]
    moc_offsets_factor = [0.0, 4.0, 8.0, 12.0]
    levels_elicitor = 20.0:5.0:100.0

    # Create figure and axes
    fig = Figure(; size=(200*length(moc_weights), 200))
    _, axs = axis_grid(
        moc_weights, 
        fig; 
        formatter=x -> "Weight = $x", 
        xlabel="Elicitor level (dB SPL)",
        ylabel="ΔL (dB)"
    )

    # Loop over parameters and generate
    lns = map(zip(axs, moc_weights)) do (ax, weight)
        # Add supplemental lines
        ablines!(ax, [-20.0], [1.0]; color=:lightgray, linestyle=:dash)
        ablines!(ax, [-30.0, -40.0, -50.0, -60.0], [1.0, 1.0, 1.0, 1.0]; color=:gray, linestyle=:dot, linewidth=0.5)

        # Plot reference line for largest expected RLF shift
        model1 = fetch_model_baseline(cf; cohc=1.0)
        sim1 = RLF(model1)
        model2 = fetch_model_baseline(cf; cohc=0.01)
        sim2 = RLF(model2)
        hlines!(ax, [threshold_shift(sim1, sim2; config=config)]; color=:gray)

        # Loop over offsets and plot each curve
        map(zip(moc_offsets_factor, colorschemes[:Dark2_8])) do (offset, c)
            # Run sim
            curve = map(levels_elicitor) do level_elicitor
                model = fetch_model(cf; moc_weight=weight, moc_offset=weight*offset)
                sim1 = RLF(model, ProbePureToneElicitor; level_elicitor=-Inf, level_spacing=1.0)
                sim2 = RLF(model, ProbePureToneElicitor; level_elicitor=level_elicitor, levels=level(sim1))
                threshold_shift(sim1, sim2; config=config)
            end

            # Plot results
            scatter!(ax, levels_elicitor, curve; color=c)
            lines!(ax, smooth(levels_elicitor, curve, :loess)...; color=c)
        end
    end

    # Add legend
    Legend(fig[1, length(moc_weights)+1], lns[1], string.(moc_offsets_factor), "Offset")

    # Adjust limits
    ylims!.(axs, 0.0, 50.0)

    # Display and save
    display(fig)
end

"""
    genfig_supp_cas_ΔL_vs_level_β([; cf=2e3, config=Config()])

Plot ΔL vs elicitor level curves for various values of MOC nonlinearity beta parameter

The goal of this figure is to elucidate how/whether β is separable from the weight and 
offset parameters.
"""
function genfig_supp_cas_ΔL_vs_level_β(
    moc_betas = [0.01, 0.02, 0.03],
    moc_weights=[10.0, 20.0];
    cf=2e3,
    config=Config(),
    plot_fit=true,
    level_step_elicitor=5.0,
    level_step_probe=2.5,
    sz_h=300,
    sz_v=150,
)  
    # Set up axes for weight sweep
    levels_elicitor = 20.0:level_step_elicitor:100.0

    # Create figure and axes
    params, fig, axs = axis_grid(
        moc_weights,
        moc_betas;
        formatter_col=x -> "β = $x", 
        formatter_row=x -> "weight = $x", 
        xlabel="Elicitor level (dB SPL)",
        ylabel="ΔL (dB)",
        sz_h=sz_h,
        sz_v=sz_v,
    )

    # Loop over parameters and generate
    map(zip(axs, params)) do (ax, (weight, β))
        # Add supplemental lines
        ablines!(ax, -10.0 .- 10.0 .* (1:1:10), ones(10); color=:gray, linestyle=:dot, linewidth=0.5)

        # Plot reference line for largest expected RLF shift
        model1 = fetch_model_baseline(cf; cohc=1.0)
        sim1 = RLF(model1)
        model2 = fetch_model_baseline(cf; cohc=0.0)
        sim2 = RLF(model2)
        hlines!(ax, [threshold_shift(sim1, sim2; config=config)]; color=:gray)

        # Loop over levels and plot curve
        curve = map(levels_elicitor) do level_elicitor
            model = fetch_model(cf; moc_weight=weight, moc_offset=weight*4.0, moc_beta=β)
            sim1 = RLF(model, ProbePureToneElicitor; level_elicitor=-Inf, level_spacing=level_step_probe)
            sim2 = RLF(model, ProbePureToneElicitor; level_elicitor=level_elicitor, levels=level(sim1))
            threshold_shift(sim1, sim2; config=config)
        end

        # Plot results
        scatter!(ax, levels_elicitor, curve; color=:black)
        if plot_fit
            lines!(ax, smooth(levels_elicitor, curve, [-10.0, 1/5, 30.0], :logistic)...; color=:black)
        end
    end

    # Adjust limits
    ylims!.(axs, 0.0, 50.0)

    # Display and save
    display(fig)
end


function genfig_survey_moc_nonlin(; config=Config())
    # Set params
    moc_weights = 5.0:5.0:15.0
    moc_offsets = 0.0:8.0:20.0
    moc_betas = 0.01 .* 2.0 .^ (-1:1:2) 

    # Create fig and axes
    params, fig, axs = axis_grid(
        moc_weights, 
        moc_offsets; 
        sz_h=350, 
        sz_v=300, 
        xlabel="MOC rate (sp/s)", 
        ylabel="COHC",
        formatter_col=x -> "Offset = $x",
        formatter_row=x -> "Weight = $x",
    )

    # Run simulation to create reference color axis for MOC rate
    model = fetch_model_baseline(2e3; stage="lsr")
    levels = 0.0:5.0:90.0
    r = @memo config simulate(RLF(model, ProbePureToneElicitor; levels=levels))

    # Run simulation to create reference color axis for COHC
    cohcs = 0.0:0.1:1.0
    ΔL = map(cohcs) do cohc
        model1 = fetch_model_baseline(2e3; cohc=1.0)
        model2 = fetch_model_baseline(2e3; cohc=cohc)
        midpoint_shift(RLF(model1, ProbePureToneElicitor), RLF(model2, ProbePureToneElicitor); config=config)
    end

    map(zip(params, axs)) do ((weight, offset), ax)
        # Plot MOC io nonlinearities 
        x = LinRange(0.0, 80.0, 100)
        map(zip(moc_betas, get.(Ref(ColorSchemes.Purples), LinRange(0.2, 1.0, length(moc_betas))))) do (β, c)
            lines!(ax, x, moc_nonlin.(x; w=weight, θ=offset*weight, β=β, minval=0.0); color=c)
        end
        
        # Plot reference data to interpret MOC rate in terms of level
        scatter!(ax, r, 1.1 .* ones(length(r)); marker=:rect, color=get.(Ref(ColorSchemes.batlow), levels ./ maximum(levels)))

        # Plot reference data to interpret COHC in terms of ΔL
        scatter!(ax, 90.0 .* ones(length(ΔL)), cohcs; marker=:rect, color=get.(Ref(reverse(ColorSchemes.glasgow)), cohcs))
    end

    # Add colorbar
    Colorbar(fig[1, length(moc_offsets)+2], limits=extrema(levels), colormap=ColorSchemes.batlow, label="Probe level (dB SPL)")
    Colorbar(fig[2, length(moc_offsets)+2], limits=extrema(ΔL), colormap=ColorSchemes.glasgow, label="ΔL (dB)")
    Colorbar(fig[3, length(moc_offsets)+2], limits=extrema(moc_betas), colormap=ColorSchemes.Purples, label="β")

    # Adjust limits
    ylims!.(axs, 0.0, 1.13)
    fig
end


"""
"""
function genfig_supp_elicitor_vs_probe_response_dynamics(; config=Config(), β=0.015)
    # Select model
    model = fetch_model(8e3; moc_beta=β)

    # Determine ipsi RLF midpoint and CAS threshold
    # These values will be used to determine the elicitor and probe levels.
    mp = midpoint(RLF(ipsi(model)); config=config)
    s_th = threshold(SLF(model; level_probe=mp); config=config)

    # Pre-run responses
    incs = -30.0:10.0:40.0
    resps = map(incs) do inc
        stim = ProbePureToneElicitor([8e3, mp], [8e3, s_th + inc]; fs=100e3)
        resp_contra = @memo config Utilities._compute(contra(model), contra(stim))
        resp_ipsi = @memo config compute(model, stim)
        y = resp_ipsi
        x = resp_contra["lsr"][contra(model).coi[1]][1:length(y)]
        μ_x = mean(x[idxswin(stim)])
        μ_y = mean(y[idxswin(stim)])
        x = x[1:100:end]
        y = y[1:100:end]
        (x, y, μ_x, μ_y)
    end

    # Create figure and axes
    fig = Figure(; size=(500, 400))
    ax = Axis(fig[1, 1])

    # Select colors
    colors = LinRange(0.0, 360.0-360.0/(length(incs)+1), length(incs))

    # Highlight maxima
    maxima = map(resps) do resp
        idx = argmax(resp[2])
        (resp[1][idx], resp[2][idx])
    end
    lines!(ax, getindex.(maxima, 1), getindex.(maxima, 2); color=:gray, linestyle=:dash)
    lines!(ax, getindex.(resps, 3), getindex.(resps, 4); color=:gray, linestyle=:dash)

    # Plot results
    s = map(enumerate(zip(resps, colors))) do (idx, (resp, c))
        # Generate colormap
        x, y, μ_x, μ_y = resp
        cm = sequential_palette(c, length(x)) 

        # Plot lines underneath
        lines!(ax, x, y; color=cm)

        # Put square marker at onset to elicitor
        idx_onset = argmax(x)
        s = scatter!(ax, [x[idx_onset]], [y[idx_onset]]; color=cm[idx_onset], markersize=15.0, marker=:rect)

        # Put circle marker at onset to probe
        idx_onset = argmax(y)
        scatter!(ax, [x[idx_onset]], [y[idx_onset]]; color=cm[idx_onset], markersize=15.0, marker=:circle)

        ax.xlabel = "Firing rate — contra elicitor LSR (sp/s)"
        ax.ylabel = "Firing rate — ipsi probe HSR (sp/s)"
        # Colorbar(fig[1, idx+1]; limits=(0.0, 0.2), colormap=cm, ticklabelsvisible=idx < length(incs) ? false : true)

        return s
    end

    # Add supporting text
    text!(ax, [55.0], [800.0]; text="Probe level = $mp dB SPL", color=:black, fontsize=12.0)
    text!(ax, [55.0], [750.0]; text="CAS threshold = $s_th dB SPL", color=:black, fontsize=12.0)

    Legend(fig[1, 2], s, string.(incs .+ s_th))
    xlims!(ax, 0.0, 120.0)
    ylims!(ax, 0.0, 900.0)
    fig
end

"""
    genfig_supp_cas_suppression_growth([; cf=2e3, config=Config()])

Plot ΔL vs elicitor level.
"""
function genfig_supp_cas_suppression_growth(
    params::MOCParams3=fetch_params("default");
    cf=2e3,
    config=Config(),
    fig=Figure(; size=(400, 300)),
    ax=Axis(fig[1, 1]),
    level_step=2.5,
    color=:black,
)
    # Pick elicitor levels
    levels_elicitor = 20.0:5.0:100.0

    # Add supplemental lines
    ints = 0.0:-10.0:-60.0
    ablines!(ax, ints, ones(length(ints)); color=:gray, linestyle=:dot, linewidth=0.5)

    # Plot reference line for largest expected RLF shift from model
    model1 = fetch_model_baseline(cf; cohc=1.0)
    sim1 = RLF(model1)
    model2 = fetch_model_baseline(cf; cohc=0.01)
    sim2 = RLF(model2)
    ΔL_max = threshold_shift(sim1, sim2; config=config)
    hlines!(ax, [threshold_shift(sim1, sim2; config=config)]; color=:gray)
    text!(ax, [3.0], [ΔL_max]; text="Maximum ΔL (model)", color=:gray)

    # Plot reference line for largest expected RLF shift from physio data
    arrows!(ax, [0.0], [12.0], [85.0], [0.0]; color=:pink)
    text!(ax, [3.0], [12.0]; text="Maximum reported ΔL (physio)", color=:pink)

    # Plot reference line for expected CAS threshold
    model = fetch_model(params, cf)
    level_probe = midpoint(RLF(ipsi(model); levels=0.0:level_step:60.0); config=config)
    slf = SLF(model; level_elicitor_step=level_step, level_probe=level_probe)
    vlines!(ax, [threshold(slf; config=config)]; color=:gray)

    # Run sim and plot
    curve = @showprogress map(levels_elicitor) do level_elicitor
        sim1 = RLF(model, ProbePureToneElicitor; level_elicitor=-Inf, level_step=level_step)
        sim2 = RLF(model, ProbePureToneElicitor; level_elicitor=level_elicitor, levels=level(sim1))
        threshold_shift(sim1, sim2; config=config)
    end
    scatter!(ax, levels_elicitor, curve; color=color)
    lines!(ax, smooth(levels_elicitor, curve)...; color=color)

    # Adjust limits
    ylims!(ax, 0.0, 50.0)
    xlims!(ax, 0.0, 100.0)
    ax.xticks = 0.0:20.0:100.0

    # Add labels
    ax.xlabel = "Elicitor level (dB SPL)"
    ax.ylabel = "ΔL (dB)"

    # Display and return
    display(fig)
    fig, ax
end

# function genfig_supp_cas_suppression_growth(
#     weight::Float64=5.0,
#     offset::Float64=5.0,
#     β::Float64=0.015;
#     kwargs...
# )
#     # Pack into MOCParams struct and then pass to method
#     params = MOCParams3(; moc_weight_func=x->weight, moc_offset=offset, moc_beta=β)
#     genfig_supp_cas_suppression_growth(params; kwargs...)
# end


"""
    genfig_supp_level_growth_panel([; config=Config()])

Survey figure depicting change in response over level in several pre-MOC stages vs CF
"""
function genfig_supp_level_growth_panel(;
    config=Config(),
    species="cat",
    cfs = round.(LogRange(0.5e3, 16e3, 41)),
    plot_refs=false,
    fig = Figure(; size=(1000, 240)),
    xlims=(0.5, 23.0),
    highlight_inset=true,
    fn="fig_supplemental_growth_curves.png",
    xticks=[],
    xminorticksvisible=true,
    ylabels=["RMS (a.u.)", "RMS (a.u.)", "Firing rate (sp/s)", "Firing rate (sp/s)"],
)  
    # Create figure object
    stages = ["c1", "ihc", "hsr", "lsr"]
    metrics = [rms, rms, mean, mean]
    axs = [tradlogxax(fig[1, i]; yscale=(i <=2 ? log10 : identity), xminorticksvisible=xminorticksvisible) for i in 1:length(stages)]
    if !isempty(xticks)
        [ax.xticks = xticks for ax in axs]
    end

    # Pick which levels will get labeled in which stages
    labelinfo = Dict(
        "c1" => [0.0, 20.0, 40.0, 80.0],
        "ihc" => [0.0, 20.0, 40.0, 80.0],
        "hsr" => [0.0, 10.0, 20.0, 30.0],
        "lsr" => [20.0, 40.0, 60.0, 80.0],
    )

    # Select levels and colors
    lvls = -15.0:5.0:80.0
    colors = get.(Ref(colorschemes[:viridis]), LinRange(0.0, 1.0, length(lvls)))

    # Map over CFs and simulate all rates
    for (idx_stage, (stage, metric, ylabel)) in enumerate(zip(stages, metrics, ylabels))
        # Select axis
        ax = axs[idx_stage]

        # Add box to indicate inset region
        if highlight_inset
            vspan!(ax, [3.0], [6.0]; color=:lightgray)
        end

        # Grab RLFs, using stage unless stage is "gain" in which case we compute gain using
        # MOC input output nonlinearity
        rlfs = map(cfs) do cf
            # Compile RLF
            rlf = RLF(
                fetch_model_baseline(cf; species=species, stage=(stage == "gain") || (stage == "gaindb") ? "mocwdr" : stage), 
                ProbePureToneElicitor; 
                levels=lvls, 
                level_elicitor=-Inf, 
                summaryfunc=metric,
            )
            @memo config simulate(rlf)
        end

        # Loop over levels and plot
        for (idx_level, (level, color)) in enumerate(zip(lvls, colors))
            # Plot iso-level curves
            scatter!(ax, cfs ./ 1e3, getindex.(rlfs, idx_level); color=color, markersize=4.0)
            if stage == "gaindb"
                x = cfs ./ 1e3
                y = getindex.(rlfs, idx_level)
                x̂, ŷ = smooth(log.(x), y, :loess)
                lines!(ax, exp.(x̂), ŷ; color=color)
            else
                lines!(ax, cfs ./ 1e3, getindex.(rlfs, idx_level); color=color)
            end

            # Plot text labels to help out
            if in(level, labelinfo[stage])
                text!(
                    ax, 
                    [1.1*cfs[end]/1e3], 
                    [rlfs[end][idx_level]]; 
                    text="$(Int(round(level)))", 
                    color=color,
                    align=(:left, :center),
                )
            end
        end
        ax.title = uppercase(stage)
        ax.ylabel = ylabel
    end
    xlims!.(axs, xlims...)

    # If plot_refs, draw vlines at 0.6 kHz
    if plot_refs
        map(axs) do ax
            vlines!(ax, [0.6]; color=:gray, linestyle=:dash)
        end
    end

    # Add labels
    Label(fig[2, :], "CF (kHz)"; fontsize=12.0)
    rowgap!(fig.layout, 1, Relative(0.0))

    # Add colorbar
    Colorbar(fig[1, 5]; limits=extrema(lvls), colormap=:viridis, label="Sound level (dB SPL)")
    
    # Display and save
    display(fig)
    saveplot(fn, fig)

    fig
end

"""
    genfig_supp_rlf_steepness([; config=Config()])

Survey figure depicting change in RLF shape across CF
"""
function genfig_supp_rlf_steepness(;
    config=Config(),
    species="cat",
    cfs = round.(LogRange(0.5e3, 16e3, 7)),
    colorscheme=:batlow,
)  
    # Create figure and axes
    fig = Figure(; size=(300, 200))
    ax = Axis(fig[1, 1])

    # Select levels and colors
    lvls = -15.0:2.0:70.0
    colors = get.(Ref(colorschemes[colorscheme]), LinRange(0.0, 1.0, length(cfs)))

    # Grab RLFs, using stage unless stage is "gain" in which case we compute gain using
    # MOC input output nonlinearity
    map(zip(cfs, colors)) do (cf, color)
        # Fetch RLF
        rlf = RLF(
            fetch_model_baseline(cf; species=species, stage="hsr");
            levels=lvls, 
        )
        r = @memo config simulate(rlf)

        # Extract midpoint, minimum, maximum
        mp = midpoint(rlf; config=config)
        minval = minimum(r)
        maxval = maximum(r)

        # Plot centered on midpoint and normalize in terms of driven rate
        lines!(ax, lvls .- mp, (r .- minval) ./ (maxval - minval); color=color)
    end

    # Add colorbar
    # TODO: use driven
    Colorbar(fig[1, 2]; limits=extrema(cfs ./ 1e3), colormap=colorscheme, label="CF (kHz)", scale=log, ticks=[1, 2, 4, 8], minorticks=IntervalsBetween(9))
    
    # Display and save
    xlims!(ax, -25.0, 25.0)
    ax.xticks = -30.0:10.0:30.0
    ax.ylabel = "Driven rate (a.u.)"
    ax.xlabel = "Level re: midpoint (dB)"
    display(fig)
    saveplot("fig_supplemental_rlf_steepness.png", fig)

    fig
end



"""
    genfig_supp_ΔL_vs_threshold_params([; cf=2e3, config=Config()])

Shows ΔL vs threshold for different parameters provided by `fetch_params`.
"""
function genfig_supp_ΔL_vs_threshold_params(
    cf=2e3; 
    config=Config(),
    level_step_probe=2.5,
    level_step_elicitor=2.5,
    level_elicitor=72.5,
    fig=Figure(; size=(2.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
)
    # Gather all params
    params = fetch_params.(["default", "current", "low weight high slope", "high weight low slope", "high weight low slope with offset", "ultrashallow"])

    # Loop over params and run at high resolution to obtain 
    results = map(params) do param
        # Fetch model
        model = fetch_model(param, cf)

        # Determine CAS threshold
        mp = midpoint(RLF(ipsi(model); level_step=level_step_probe); config=config)
        cas_threshold = threshold(
            SLF(model; level_probe=mp, level_elicitor_step=level_step_elicitor); 
            config=config,
        )

        # Determine ΔL at fixed elicitor level
        cas_ΔL = threshold_shift(
            RLFPair(model, [-Inf, level_elicitor]; level_step=level_step_probe)...; 
            config=config
        )
        (cas_threshold, cas_ΔL)
    end

    # Fetch all thresholds and ΔLs
    ths = getindex.(results, 1)
    ΔLs = getindex.(results, 2)

    # Create figure and plot
    scatter!(ax, ths, ΔLs)

    # Labels and lims
    ax.xlabel = "CAS threshold (dB SPL)"
    ax.ylabel = "ΔL (dB)"
    ylims!(ax, 0.0, 20.0)
    xlims!(ax, 30.0, 60.0)

    # Display and return
    display(fig)
    fig
end


"""
    genfig_supp_ΔL_vs_threshold_grid([; cf=2e3, config=Config()])

Shows ΔL vs threshold for different parameters provided by `fetch_params`.
"""
function genfig_supp_ΔL_vs_threshold_grid(
    weights=LinRange(1.0, 40.0, 9),
    offsets=LinRange(0.0, 40.0, 9),
    betas=LogRange(0.001, 0.05, 9);
    cf=2e3, 
    config=Config(),
    level_step_probe=1.0,
    level_step_elicitor=1.0,
    level_elicitor=72.5,
    fig=Figure(; size=(4.5*100, 2*100)),
    ax=Axis(fig[1, 1]),
    ax_supp=Axis(fig[1, 2]),
    colors=get.(Ref(colorschemes[:viridis]), LinRange(0.0, 1.0, length(betas))),
)
    # Annotate box to indicate rough position of expected data
    lb_th, ub_th = 
        query_cas_threshold_env_Warren1989a_Fig6(cf),
        query_cas_threshold_Warren1989a_Fig6(cf)
    vspan!(ax, [lb_th], [ub_th]; color=(:gray, 0.2))

    lb_ΔL, ub_ΔL =
        query_cas_ΔL_Warren1989b_Fig2(cf),
        query_cas_ΔL_env_Warren1989b_Fig2(cf)
    hspan!(ax, [lb_ΔL], [ub_ΔL]; color=(:gray, 0.2))

    # Gather all params
    params = map(Iterators.product(weights, offsets, betas)) do (weight, offset, β)
        MOCParams3(
            moc_weight_func=cf -> cf < 2e3 ? weight : weight * peaknorm_gaussian(log2(cf/2e3), 0.0, 2.0),
            moc_offset=offset, 
            moc_beta=β
        )
    end

    # Loop over params and run at high resolution to obtain 
    results = @showprogress map(params) do param
        # Fetch model
        model = fetch_model(param, cf)

        # Determine CAS threshold
        mp = midpoint(RLF(ipsi(model); level_step=level_step_probe); config=config)
        cas_threshold = threshold(
            SLF(model; level_probe=mp, level_elicitor_step=level_step_elicitor); 
            config=config,
        )

        # Determine ΔL at fixed elicitor level
        cas_ΔL = threshold_shift(
            RLFPair(model, [-Inf, level_elicitor]; level_step=level_step_probe)...; 
            config=config
        )
        (cas_threshold, cas_ΔL)
    end

    # Flatten results into vector
    params = vcat(params...)
    color_idxs = map(getfield.(params, :moc_beta)) do β
        findall(β .== betas)[1]
    end
    results = vcat(results...)

    # Fetch all thresholds and ΔLs
    ths = getindex.(results, 1)
    ΔLs = getindex.(results, 2)

    # Create figure and plot
    scatter!(ax, ths, ΔLs; color=colors[color_idxs])

    # Labels and lims
    ax.xlabel = "CAS threshold (dB SPL)"
    ax.ylabel = "ΔL (dB)"
    ylims!(ax, 0.0, 50.0)
    xlims!(ax, 10.0, 80.0)

    # Print to console all params that satisfy boundaries
    idxs = findall((lb_th .< ths) .& (ths .< ub_th) .& (lb_ΔL .< ΔLs) .& (ΔLs .< ub_ΔL))
    println("Found $(length(idxs)) parameters that satisfy boundaries")
    for idx in idxs
        println(params[idx])
    end

    # Display weight at CF and β that yield good results
    scatter!(
        ax_supp, 
        getfield.(params[idxs], :moc_beta),
        map(x -> x(cf), getfield.(params[idxs], :moc_weight_func));
    )
    ax_supp.xlabel = "Beta"
    ax_supp.ylabel = "Weight"
    ax_supp.xscale = log10

    # Determine "best" fit as fit closest to center of ΔL
    dist = sqrt.((ΔLs .- mean([lb_ΔL, ub_ΔL])) .^2  + (ths .- mean([lb_th, ub_th])) .^2)
    idxs_valid = .!isnan.(dist)
    best = params[idxs_valid][argmin(dist[idxs_valid])]

    # Display and return
    display(fig)
    params[idxs], best
end

function genfig_supp_tc_diagnose(; 
    config=Config(),
    size=(1.6*100, 1.3*100),
    fig=Figure(; size=size),
    ax=Axis(fig[1, 1]; xscale=log10),
    level_high=100.0,
)
    # Plot
    model = fetch_model(2e3)
    sim = RLFTC(ipsi(model); level_high=level_high, freq_low=-2.0, freq_high=2.0, n_freq=61)
    x̂, ŷ = smooth(log2.(Utilities.freq(sim)), threshold_curve(sim; config=config), :loess)
    lines!(ax, (2 .^ x̂) ./ 1e3, ŷ; color=:black, linewidth=1.0)

    # Scale
    ylims!(ax, -20.0, level_high)
    xlims!(ax, 2/8, 2*8)
    ax.xticks = [1, 2, 4]

    # Add labels
    ax.xlabel = "Frequency (kHz)"
    ax.ylabel = "Threshold (dB SPL)"

    # Display and save
    display(fig)
    fig
end


"""
    genfig_supp_c1_compression([; config=Config()])
"""
function genfig_supp_c1_compression(;
    config=Config(),
    species="cat",
    cfs = round.(LogRange(0.25e3, 16e3, 10)),
)  
    # Select levels and colors
    lvls = -15.0:5.0:80.0
    colors = get.(Ref(colorschemes[:viridis]), LinRange(0.0, 1.0, length(cfs)))

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1]; yscale=log10)

    # Map over CFs, simulate response, and plot
    lns = map(zip(cfs, colors)) do (cf, c)
        # Compile RLF
        rlf = RLF(
            fetch_model_baseline(cf; species=species, stage="c1"),
            ProbePureToneElicitor; 
            levels=lvls, 
            level_elicitor=-Inf, 
            summaryfunc=rms,
        )
        curve = @memo config simulate(rlf)

        scatter!(ax, lvls, curve .+ eps(); color=c)
        lines!(ax, smooth(lvls, curve .+ eps(), :loess)...; color=c)
    end
    Legend(fig[1, 2], lns, string.(round.(cfs ./ 1e3; digits=2)); title="CF (kHz)")

    ax.xlabel = "Sound level (dB SPL)"
    ax.ylabel = "C1 response (rms)"

    fig
end


function genfig_supp_RLFs_debug(; config=Config(), fiber_type="hsr")
    # Choose CFs
    cfs = LogRange(0.5e3, 32e3, 25)

    # Create figure and axes
    fig = Figure(; size=(1000, 1000))
    n_ax_side = Int(ceil(sqrt(length(cfs))))
    axs = [Axis(fig[i, j]) for i in 1:n_ax_side, j in 1:n_ax_side]

    # Fetch models and simulations
    models = map(x -> fetch_model_baseline(x; stage=fiber_type), cfs)
    if fiber_type == "hsr"
        lvls = -5.0:1.0:45.0 
    else
        lvls = 10.0:1.0:80.0
    end
    sims = RLF.(models; levels=lvls)
    colors = get.(Ref(colorschemes[:viridis]), LinRange(0.0, 1.0, length(cfs)))

    # Loop over CFs, simulate response, and plot
    for (sim, color, ax) in zip(sims, colors, axs)
        # Fetch and response
        r = @memo config simulate(sim)

        # Extract MP
        mp = Utilities.nthpoint(level(sim), r, 0.5)

        # Extract threshold
        th = Utilities.threshold(level(sim), r)

        # Plot RLF
        lines!(ax, level(sim), r; color=color, linewidth=2.0)

        # Plot midpoint
        scatter!(ax, [th], [r[argmin(abs.(level(sim) .- th))]]; color=color, markersize=10.0)
        scatter!(ax, [mp], [r[argmin(abs.(level(sim) .- mp))]]; color=color, markersize=10.0)
        vlines!(ax, [mp]; color=color, linestyle=:dash, linewidth=0.5)

        # Set limits and title
        ylims!(ax, 0.0, fiber_type == "hsr" ? 450.0 : 150.0)
        ax.title = string("CF: ", round(getcf(sim) / 1e3, sigdigits=2), " kHz")
    end
    [ax.yticklabelsvisible = false for ax in axs[:, 2:end]]
    [ax.xticklabelsvisible = false for ax in axs[1:(size(axs)[2]-1), :]]
    fig
end