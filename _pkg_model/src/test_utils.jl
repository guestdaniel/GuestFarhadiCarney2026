export testparams, postprocess_simulations, run_2014_vs_2023, run_2023_vs_2023, rtol_2014_vs_2023, rtol_2023_vs_2023, plot2

# Define peripheral and post-peripheral stages to test
stages_peripheral = ["control", "c1", "c2", "ihc", "hsr", "lsr"]
stages_peripheral_multichannel = ["control", "c1", "c2", "ihc"]
stages_subcortical = ["cn", "ic"]

# Define parameter sets at which we'll test SFIE implementations
params_sfie = [
    (1.0e-3, 2e-3, 1e-3, 1.0, 0.5),  # τ_e, τ_i, d, a, s
    (1.5e-3, 1e-3, 2e-3, 2.0, 1.2),
]

# Dictionary with above data
testparams = Dict(
    "stages_peripheral" => stages_peripheral,
    "stages_subcortical" => stages_subcortical,
    "stages_peripheral_multichannel" => stages_peripheral_multichannel,
    "params_sfie" => params_sfie,
)

"""
    plot2(one, two)
"""
function plot2(one, two)
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, timevec(one, 100e3), one)
    lines!(ax, timevec(two, 100e3), two)
    ax = Axis(fig[2, 1])
    lines!(ax, timevec(one, 100e3), two .- one)
    display(fig)
    fig, ax
end

"""
    rtol_2014_vs_2023(stage)

Get suitable tolerance value for a given stage
"""
function rtol_2014_vs_2023(stage)
    @match stage begin
        "syn" => 0.08
        "hsr" => 0.08
        "lsr" => 0.08
        _ => 0.01
    end
end

"""
    rtol_2014_vs_2023(stage)

Get suitable tolerance value for a given stage
"""
function rtol_2023_vs_2023(stage)
    @match stage begin
        _ => 0.02
    end
end

"""
    postprocess_simulations(sim::Vector, stage::String, model::string, cf::Float64)

Postprocesses a simulation based on which model produced it to allow for testing

# Arguments
- `sim`: Vector containing some simulation result for a given stage
- `stage`: String indicating which stage this stage corresponds to, from ["control", "c1", \
    "c2", "ihc", "expon", "sout1", "syn", "hsr", "lsr", "cn"]
- `model`: Which model the response came from, from ["zbc2014", "gfc2023"]
- `cf`: Characteristic frequency (Hz)
"""
function postprocess_simulations(
    sim::Vector, 
    stage::String, 
    model::String, 
    cf::Float64;
    downsample=true,
    avoid_irregularities=true,
)
    # If we're looking at control, c1, or c2 in an old model, we need to shift signal by 
    # delaypoint samples to accomodate the fact that we shifted delay from immediately after
    # IHC in original code to immediately after middle ear filter in new code
    if (model == "zbc2014") & (stage in ["control", "c1", "c2"])
        delay = ccall(
            (:delay_cat, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
            Cdouble,
            (
                Cdouble,
            ),
            cf,
        )
        delaypoint = Int(round(delay/(1/100e3)))
        sim = shiftsignal(sim, delaypoint)
    end

    # If we're looking at power-law synapse stages, we need to adjust for the fact that
    # the new model operates at a single continuous sampling rate, while the old model
    # operates at a lower internal sampling rate for the power-law synapse stage 
    # To account for this, we downsample the new outputs by simply selecting every
    # 10th sample. We also need to account for the fact that the original sout1 and 
    # sout2 (and the whole powerlaw in general) included zeropadding on the edges of
    # the IHC respons. This is eliminated in the new code and produces edge effects numerical
    # the temporal edge of the simulations.
    if (model == "gfc2023") & (stage in ["sout1_hsr", "sout2_hsr", "sout1_lsr", "sout2_lsr"])
        if downsample; sim = sim[1:10:end]; end;
    end
    if (model == "zbc2014") & (stage in ["sout1", "sout2"])
        delaypoint = Int(round(7500/(cf/1e3)))
        delaypoint = Int(ceil(delaypoint/10))
        sim = sim[(1:(length(sim) - delaypoint*2)) .+ delaypoint]
    end

    # If we're looking at the synapse/rate, we need to get every 10th sample (this is because 
    # the original code used a linear interpolation to upsample back to the stimulus 
    # sampling rate, but the new code is actually simulated at 100 kHz, producing large
    # disparities between sample points)
    if stage in ["syn_hsr", "syn_lsr", "syn", "hsr", "lsr"]
        if downsample; sim = sim[1:10:end]; end;
    end

    # Return (possibly subsetted) data
    if avoid_irregularities
        if stage == "control"
            # Control onset is messed up a bit because it starts out non-zero, so simply 
            # zero-padding old control signal isn't viable and we only want to look at the 
            # relevant pieces
            sim = sim[2000:end]
        elseif stage in ["sout1_hsr", "sout2_hsr", "syn_hsr", "sout1_lsr", "sout2_lsr", "syn_lsr", "sout1", "sout2", "syn", "hsr", "lsr"]
            # For synapse stuff, we need to avoid the initial few samples because the lack of
            # a "delaypoint" system in the new model creates onset irregularities
            sim = sim[50:end]
        end
    end

    return sim
end

"""
    postprocess_simulations(sim::Dict, model::String, cf::Vector{Float64})

Postprocesses a full simulation (in format of output of sim_gfc2023_dict)
"""
function postprocess_simulations(sim::Dict, model::String, cf::Vector{Float64}; kwargs...)
    sim = map(collect(keys(sim))) do key
        key => map(zip(sim[key], cf)) do (chan, _cf)
            postprocess_simulations(chan, key, model, _cf; kwargs...)
        end
    end
    Dict(sim...)
end

"""
    run_2014_vs_2023(x, cf)

Simulates responses for old and new model at a stage for a short pure tone stimulus
"""
function run_2014_vs_2023(
    x::Vector{Float64}, 
    cf::Vector{Float64},
    args_orig=Dict{Symbol, Any}(),
    args_new=Dict{Symbol, Any}(); 
    kwargs...
)
    # Simulate 2014 response
    orig = map(f -> sim_orig_dict(x, f; args_orig...), cf)

    # Restructure 2014 response to emulate 2023 response layout (dictionary of vector of vector)
    temp = map(collect(keys(orig[1]))) do key
        key => [x[key] for x in orig]
    end
    orig = Dict(temp...)

    # Simulate 2023 response
    new = sim_gfc2023_dict(x, cf; args_new...)

    # Loop through each response type and channel and postprocess
    orig = postprocess_simulations(orig, "zbc2014", cf; kwargs...)
    new = postprocess_simulations(new, "gfc2023", cf; kwargs...)

    return orig, new
end

"""
    run_2023_vs_2023(x, cf)

Simulates responses for new model with different parameters
"""
function run_2023_vs_2023(
    x::Vector{Float64}, 
    cf::Vector{Float64},
    args_a=Dict{Symbol, Any}(),
    args_b=Dict{Symbol, Any}(); 
    kwargs...
)
    a = sim_gfc2023_dict(x, cf; args_a...)
    a = postprocess_simulations(a, "gfc2023", cf; kwargs...)
    b = sim_gfc2023_dict(x, cf; args_b...)
    b = postprocess_simulations(b, "gfc2023", cf; kwargs...)

    return a, b
end