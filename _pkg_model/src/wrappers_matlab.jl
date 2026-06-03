export sim_gfc2023_wrapper_mex, sim_gfc2023_wrapper_dict_mex

function sim_gfc2023_wrapper_mex(
    x::Vector{Float64}, 
    cf::Vector{Float64}; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    fractional=false,
    ic_tau_e=1.0e-3,
    ic_tau_i=2.0e-3,
    ic_delay=1.0e-3,
    ic_amp=4.0,
    ic_inh=0.9,
    moc_cutoff=0.2,
    moc_beta_wdr=0.01,
    moc_offset_wdr=0.0,
    moc_beta_ic=0.01,
    moc_offset_ic=0.0,
    moc_weight_wdr=0.0,
    moc_weight_ic=0.0,
    moc_width_wdr=0.5,
)
    # Calculate n_chan
    n_chan = length(cf)
    len_total = length(x)

    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    # Synthesize ffGn
    if fractional
        ffGn_hsr = map(1:n_chan) do _
            ffGn_native(
                length(x),
                1/fs,
                0.9,
                1.0,
                100.0,
            )
        end
        ffGn_lsr = map(1:n_chan) do _
            ffGn_native(
                length(x),
                1/fs,
                0.9,
                1.0,
                0.1,
            )
        end
    else
        ffGn_hsr = [zeros(len_total) for _ in 1:n_chan]
        ffGn_lsr = [zeros(len_total) for _ in 1:n_chan]
    end
    ffGn_hsr = permutedims(hcat(ffGn_hsr...))
    ffGn_lsr = permutedims(hcat(ffGn_lsr...))

    # Re-represent x as a matrix instead of vector
    x = permutedims(reshape(x, (size(x)..., 1)))

    # Call matlab using special matlab string syntax from MATLAB.jl
    mat"""
        [$ihcout, $hsrout, $lsrout, $icout, $gain] = sim_efferent_model_mex(...
            $x, ...
            $ffGn_hsr, ...
            $ffGn_lsr, ...
            $cf, ...
            1, ...
            $(1/fs), ...
            $cohc, ...
            $cihc, ...
            1, ...
            $ic_tau_e, ...
            $ic_tau_i, ...
            $ic_delay, ...
            $ic_amp, ...
            $ic_inh, ...
            $moc_cutoff, ...
            $moc_beta_wdr, ...
            $moc_offset_wdr, ...
            $moc_beta_ic, ...
            $moc_offset_ic, ...
            $moc_weight_wdr, ...
            $moc_weight_ic, ...
            $moc_width_wdr ...
        );
    """

    # Return
    new = map([ihcout, hsrout, lsrout, icout, gain]) do output
        mapslices(x->[x], output, dims=2)[:]
    end
    return new 
end

function sim_gfc2023_wrapper_mex(x::Vector{Float64}, cf::Float64; kwargs...)
    [x[1] for (idx, x) in enumerate(sim_gfc2023_wrapper_mex(x, [cf]; kwargs...))]
end

function sim_gfc2023_wrapper_dict_mex(args...; kwargs...)
    ihc, hsr, lsr, ic, gain = sim_gfc2023_wrapper_mex(args...; kwargs...)
    return Dict(
        "ihc" => ihc,
        "hsr" => hsr,
        "lsr" => lsr,
        "ic" => ic,
        "gain" => gain,
    )
end
