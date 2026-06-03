export sim_gfc2023, sim_gfc2023!, sim_gfc2023_dict, sim_gfc2023_dict!
export GFC2023_Mem, update_ffGn!, zero_state!

# ##################################################################################################
# GFC2023_Mem structure to store pre-allocated memory for sim_gfc2023!
# ##################################################################################################

"""
    GFC2023_Mem

Structure to hold pre-allocated memory for sim_gfc2023!
"""
struct GFC2023_Mem
    # Scalar arguments that determine memory size
    fs::Float64
    len_total::Int
    n_chan::Int

    # Pre-allocated memory
    ffGn_hsr::Vector{Vector{Float64}}
    ffGn_lsr::Vector{Vector{Float64}}
    controlout::Vector{Vector{Float64}}
    c1out::Vector{Vector{Float64}}
    c2out::Vector{Vector{Float64}}
    ihcout::Vector{Vector{Float64}}
    expout_hsr::Vector{Vector{Float64}}
    sout1_hsr::Vector{Vector{Float64}}
    sout2_hsr::Vector{Vector{Float64}}
    synout_hsr::Vector{Vector{Float64}}
    expout_lsr::Vector{Vector{Float64}}
    sout1_lsr::Vector{Vector{Float64}}
    sout2_lsr::Vector{Vector{Float64}}
    synout_lsr::Vector{Vector{Float64}}
    hsrout::Vector{Vector{Float64}}
    lsrout::Vector{Vector{Float64}}
    cnout::Vector{Vector{Float64}}
    icout::Vector{Vector{Float64}}
    mocwdr::Vector{Vector{Float64}}
    mocic::Vector{Vector{Float64}}
    gain::Vector{Vector{Float64}}
    gainpostmix::Vector{Vector{Float64}}
end

# Method for GFC2023_Mem in terms of len_total and n_chan
function GFC2023_Mem(len_total::Int64, n_chan::Int64, fs::Float64=100e3)
    # Start by initializing empty zerod ffGN (fractional == false)
    ffGn_hsr, ffGn_lsr = get_ffGn(len_total, false, n_chan, fs)

    # Pre-allocate memory for intermediate/output variables; all values are initialized 
    # at 0, except for gain and gainpostmix, which are initialized at 1. 
    controlout = [zeros(len_total) for _ in 1:n_chan]
    c1out = [zeros(len_total) for _ in 1:n_chan]
    c2out = [zeros(len_total) for _ in 1:n_chan]
    ihcout = [zeros(len_total) for _ in 1:n_chan]
    expout_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_hsr = [zeros(len_total) for _ in 1:n_chan]
    synout_hsr = [zeros(len_total) for _ in 1:n_chan]
    expout_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_lsr = [zeros(len_total) for _ in 1:n_chan]
    synout_lsr = [zeros(len_total) for _ in 1:n_chan]
    hsrout = [zeros(len_total) for _ in 1:n_chan]
    lsrout = [zeros(len_total) for _ in 1:n_chan]
    cnout = [zeros(len_total) for _ in 1:n_chan]
    icout = [zeros(len_total) for _ in 1:n_chan]
    mocwdr = [zeros(len_total) for _ in 1:n_chan]
    mocic = [zeros(len_total) for _ in 1:n_chan]
    gain = [ones(len_total) for _ in 1:n_chan]
    gainpostmix = [ones(len_total) for _ in 1:n_chan]

    # Wrap everything in GFC2023_Mem struct call
    GFC2023_Mem(
        fs,
        len_total,
        n_chan,
        ffGn_hsr,
        ffGn_lsr,
        controlout,
        c1out,
        c2out,
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
        gainpostmix
    )
end

# Method for GFC2023_Mem that matches function signature of sim_gfc2023!
function GFC2023_Mem(
    x::Vector{Float64},
    cfs::Vector{Float64};
    fs=100e3,
    dur_pad_left=0.02,
    dur_pad_right=0.0,
)
    len_total = length(x) + Int(floor(dur_pad_left * fs)) + Int(floor(dur_pad_right * fs))
    n_chan = length(cfs)
    return GFC2023_Mem(len_total, n_chan, fs)
end

# Method for GFC2023_Mem that fills in ffGn values when fractional == true
function update_ffGn!(mem::GFC2023_Mem)
    # Synthesize ffGn
    for i in eachindex(mem.ffGn_hsr)
        mem.ffGn_hsr[i] .= ffGn_native(
            mem.len_total,
            1 / mem.fs,
            0.9,
            1.0,
            100.0,
        )
        mem.ffGn_lsr[i] .= ffGn_native(
            mem.len_total,
            1 / mem.fs,
            0.9,
            1.0,
            0.1,
        )
    end
end

# Method for GFC2023_Mem that zeros out state variables
function zero_state!(mem::GFC2023_Mem)
    for chan in 1:mem.n_chan
        fill!(mem.ffGn_hsr[chan], 0.0)
        fill!(mem.ffGn_lsr[chan], 0.0)
        fill!(mem.controlout[chan], 0.0)
        fill!(mem.c1out[chan], 0.0)
        fill!(mem.c2out[chan], 0.0)
        fill!(mem.ihcout[chan], 0.0)
        fill!(mem.expout_hsr[chan], 0.0)
        fill!(mem.sout1_hsr[chan], 0.0)
        fill!(mem.sout2_hsr[chan], 0.0)
        fill!(mem.synout_hsr[chan], 0.0)
        fill!(mem.expout_lsr[chan], 0.0)
        fill!(mem.sout1_lsr[chan], 0.0)
        fill!(mem.sout2_lsr[chan], 0.0)
        fill!(mem.synout_lsr[chan], 0.0)
        fill!(mem.hsrout[chan], 0.0)
        fill!(mem.lsrout[chan], 0.0)
        fill!(mem.cnout[chan], 0.0)
        fill!(mem.icout[chan], 0.0)
        fill!(mem.mocwdr[chan], 0.0)
        fill!(mem.mocic[chan], 0.0)
        fill!(mem.gain[chan], 1.0)
        fill!(mem.gainpostmix[chan], 1.0)
    end
end

# ##################################################################################################
# Fractional Gaussian noise help code
# ##################################################################################################

"""
    get_ffGn(len_total, fractional, n_chan; fs=100e3)

Function to prepare ffGn inputs for sim_gfc2023 and sim_gfc2023!
"""
function get_ffGn(len_total, fractional, n_chan, fs=100e3)
    # Synthesize ffGn
    if fractional
        ffGn_hsr = map(1:n_chan) do _
            ffGn_native(
                len_total,
                1 / fs,
                0.9,
                1.0,
                100.0,
            )
        end
        ffGn_lsr = map(1:n_chan) do _
            ffGn_native(
                len_total,
                1 / fs,
                0.9,
                1.0,
                0.1,
            )
        end
    else
        ffGn_hsr = [zeros(len_total) for _ in 1:n_chan]
        ffGn_lsr = [zeros(len_total) for _ in 1:n_chan]
    end

    return ffGn_hsr, ffGn_lsr
end

# ##################################################################################################
# Primary model wrapper code
# ##################################################################################################

""" 
    sim_gfc2023(input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates full model output for sound-pressure input

# Positional arguments 
- `x::Vector{Float64}`: sound-pressure waveform (Pa)
- `cf::Float64`: characteristic frequency of the fiber in Hz

# Keyword arguments
- `fs::Float64`: sampling rate of the *input* in Hz
- `cohc::Float64`:
- `cihc::Float64`:
- `species::String`:
- `fractional::Bool`: 

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_gfc2023(
    x::Vector{Float64},
    cf::Vector{Float64};
    gain::Vector{Vector{Float64}}=[Float64[]],
    fs::Float64=100e3,
    cohc::Vector{Float64}=ones(size(cf)),
    cihc::Vector{Float64}=ones(size(cf)),
    species::String="human",
    fractional=false,
    powerlaw_mode=2,
    cn_tau_e=0.5e-3,
    cn_tau_i=2.0e-3,
    cn_delay=1.0e-3,
    cn_amp=1.5,
    cn_inh=0.6,
    ic_tau_e=1.0 / (10.0 * 64.0),  # BMF = 64 Hz
    ic_tau_i=1.0 / (10.0 * 64.0) * 1.5,
    ic_delay=1.0 / (10.0 * 64.0) * 2.0,
    ic_amp=1.0,
    ic_inh=0.9,
    moc_cutoff=0.64,
    moc_beta=fill(0.2, length(cf)),
    moc_offset=fill(5.0, length(cf)),
    moc_minval=0.1,
    moc_maxval=1.0,
    moc_weight=fill(1.0, length(cf)),
    moc_width=0.5,
    dur_pad_left=0.02,
    moc_delay=0.025,
    moc_fix_gain=false,
    clip_left=dur_pad_left == 0.0 ? false : true,
    dur_pad_right=0.0,
    clip_right=dur_pad_right == 0.0 ? false : true,
)::Vector{Vector{Vector{Float64}}}
    # Calculate pad sizes in samples
    len_pad_left = Int(floor(dur_pad_left * fs))
    len_pad_right = Int(floor(dur_pad_right * fs))
    len_stim = length(x)
    len_total = len_pad_left + len_stim + len_pad_right

    # Pad x
    stim = vcat(zeros(len_pad_left), x, zeros(len_pad_right))

    # Calculate n_chan
    n_chan = length(cf)

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
                len_total,
                1 / fs,
                0.9,
                1.0,
                100.0,
            )
        end
        ffGn_lsr = map(1:n_chan) do _
            ffGn_native(
                len_total,
                1 / fs,
                0.9,
                1.0,
                0.1,
            )
        end
    else
        ffGn_hsr = [zeros(len_total) for _ in 1:n_chan]
        ffGn_lsr = [zeros(len_total) for _ in 1:n_chan]
    end

    # If MOC weight is passed as a scalar, replace it with a vector of the same length as cf
    # filling in the scalar weight. Same applies to moc_beta and moc_offset.
    if typeof(moc_weight) == Float64
        moc_weight = fill(moc_weight, length(cf))
    end
    if typeof(moc_beta) == Float64
        moc_beta = fill(moc_beta, length(cf))
    end
    if typeof(moc_offset) == Float64
        moc_offset = fill(moc_offset, length(cf))
    end

    # Pre-allocate memory
    controlout = [zeros(len_total) for _ in 1:n_chan]
    c1out = [zeros(len_total) for _ in 1:n_chan]
    c2out = [zeros(len_total) for _ in 1:n_chan]
    ihcout = [zeros(len_total) for _ in 1:n_chan]
    expout_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_hsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_hsr = [zeros(len_total) for _ in 1:n_chan]
    synout_hsr = [zeros(len_total) for _ in 1:n_chan]
    expout_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout1_lsr = [zeros(len_total) for _ in 1:n_chan]
    sout2_lsr = [zeros(len_total) for _ in 1:n_chan]
    synout_lsr = [zeros(len_total) for _ in 1:n_chan]
    hsrout = [zeros(len_total) for _ in 1:n_chan]
    lsrout = [zeros(len_total) for _ in 1:n_chan]
    cnout = [zeros(len_total) for _ in 1:n_chan]
    icout = [zeros(len_total) for _ in 1:n_chan]
    mocwdr = [zeros(len_total) for _ in 1:n_chan]
    mocic = [zeros(len_total) for _ in 1:n_chan]
    if isempty(gain[1])
        gain = [ones(len_total) for _ in 1:n_chan]
    end
    gainpostmix = [ones(len_total) for _ in 1:n_chan]

    # Add length assertion
    @assert length(gain) == length(cf)

    # Run model
    model!(
        stim,
        ffGn_hsr,
        ffGn_lsr,
        cf,
        n_chan,
        1 / fs,
        len_total,
        cohc,
        cihc,
        species_flag,
        100.0,
        powerlaw_mode,
        cn_tau_e,
        cn_tau_i,
        cn_delay,
        cn_amp,
        cn_inh,
        ic_tau_e,
        ic_tau_i,
        ic_delay,
        ic_amp,
        ic_inh,
        moc_cutoff,
        moc_beta,
        moc_offset,
        moc_minval,
        moc_maxval,
        moc_weight,
        moc_width,
        dur_pad_left,
        moc_delay,
        Int(moc_fix_gain),
        controlout,
        c1out,
        c2out,
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
        gainpostmix,
    )

    # Return
    outputs = [
        controlout,
        c1out,
        c2out,
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
        gainpostmix
    ]
    if clip_left | clip_right
        outputs = map(outputs) do output
            output = map(output) do channel
                idx_left = clip_left ? (len_pad_left + 1) : 1
                idx_right = clip_right ? length(channel) - len_pad_right : length(channel)
                channel = channel[idx_left:idx_right]
            end
        end
    end
    return outputs
end

""" 
    sim_gfc2023!(mem..., input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)
    sim_gfc2023!(mem::GFC2023_Mem, input, cf; fs=100e3, fs_synapse=10e3, power_law="approximate", fractional=false, n_rep=1)

Simulates full model output for sound-pressure input in place using pre-allocated memory mem

# Arguments
- `ffGn_hsr::Vector{Vector{Float64}}`:
- `ffGn_lsr::Vector{Vector{Float64}}`:
- `controlout::Vector{Vector{Float64}}`:
- `c1out::Vector{Vector{Float64}}`:
- `c2out::Vector{Vector{Float64}}`:
- `ihcout::Vector{Vector{Float64}}`:
- `expout_hsr::Vector{Vector{Float64}}`:
- `sout1_hsr::Vector{Vector{Float64}}`:
- `sout2_hsr::Vector{Vector{Float64}}`:
- `synout_hsr::Vector{Vector{Float64}}`:
- `expout_lsr::Vector{Vector{Float64}}`:
- `sout1_lsr::Vector{Vector{Float64}}`:
- `sout2_lsr::Vector{Vector{Float64}}`:
- `synout_lsr::Vector{Vector{Float64}}`:
- `hsrout::Vector{Vector{Float64}}`:
- `lsrout::Vector{Vector{Float64}}`:
- `cnout::Vector{Vector{Float64}}`:
- `icout::Vector{Vector{Float64}}`:
- `mocwdr::Vector{Vector{Float64}}`:
- `mocic::Vector{Vector{Float64}}`:
- `gain::Vector{Vector{Float64}}`:
- `gainpostmix::Vector{Vector{Float64}}`:
- `x::Vector{Float64}`: sound-pressure waveform (Pa)
- `cf::Float64`: characteristic frequency of the fiber in Hz
- `fs::Float64`: sampling rate of the *input* in Hz
- kwargs...

# Returns
- `output::Vector{Float64}`: synapse output (unknown units?), length is `length(input)`
"""
function sim_gfc2023!(
    ffGn_hsr::Vector{Vector{Float64}},
    ffGn_lsr::Vector{Vector{Float64}},
    controlout::Vector{Vector{Float64}},
    c1out::Vector{Vector{Float64}},
    c2out::Vector{Vector{Float64}},
    ihcout::Vector{Vector{Float64}},
    expout_hsr::Vector{Vector{Float64}},
    sout1_hsr::Vector{Vector{Float64}},
    sout2_hsr::Vector{Vector{Float64}},
    synout_hsr::Vector{Vector{Float64}},
    expout_lsr::Vector{Vector{Float64}},
    sout1_lsr::Vector{Vector{Float64}},
    sout2_lsr::Vector{Vector{Float64}},
    synout_lsr::Vector{Vector{Float64}},
    hsrout::Vector{Vector{Float64}},
    lsrout::Vector{Vector{Float64}},
    cnout::Vector{Vector{Float64}},
    icout::Vector{Vector{Float64}},
    mocwdr::Vector{Vector{Float64}},
    mocic::Vector{Vector{Float64}},
    gain::Vector{Vector{Float64}},
    gainpostmix::Vector{Vector{Float64}},
    x::Vector{Float64},
    cf::Vector{Float64};
    fs::Float64=100e3,
    cohc::Vector{Float64}=ones(size(cf)),
    cihc::Vector{Float64}=ones(size(cf)),
    species::String="human",
    powerlaw_mode=2,
    cn_tau_e=0.5e-3,
    cn_tau_i=2.0e-3,
    cn_delay=1.0e-3,
    cn_amp=1.5,
    cn_inh=0.6,
    ic_tau_e=1.0 / (10.0 * 64.0),  # BMF = 64 Hz
    ic_tau_i=1.0 / (10.0 * 64.0) * 1.5,
    ic_delay=1.0 / (10.0 * 64.0) * 2.0,
    ic_amp=1.0,
    ic_inh=0.9,
    moc_cutoff=0.64,
    moc_beta=fill(0.2, length(cf)),
    moc_offset=fill(5.0, length(cf)),
    moc_minval=0.1,
    moc_maxval=1.0,
    moc_weight=fill(1.0, length(cf)),
    moc_width=0.5,
    dur_pad_left=0.02,
    moc_delay=0.025,
    moc_fix_gain=false,
    clip_left=dur_pad_left == 0.0 ? false : true,
    dur_pad_right=0.0,
    clip_right=dur_pad_right == 0.0 ? false : true,
)::Vector{Vector{Vector{Float64}}}
    # Calculate pad sizes in samples
    len_pad_left = Int(floor(dur_pad_left * fs))
    len_pad_right = Int(floor(dur_pad_right * fs))
    len_stim = length(x)
    len_total = len_pad_left + len_stim + len_pad_right

    # Pad x
    stim = vcat(zeros(len_pad_left), x, zeros(len_pad_right))

    # Calculate n_chan
    n_chan = length(cf)

    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    # If MOC weight is passed as a scalar, replace it with a vector of the same length as cf
    # filling in the scalar weight. Same applies to moc_beta and moc_offset.
    if typeof(moc_weight) == Float64
        moc_weight = fill(moc_weight, length(cf))
    end
    if typeof(moc_beta) == Float64
        moc_beta = fill(moc_beta, length(cf))
    end
    if typeof(moc_offset) == Float64
        moc_offset = fill(moc_offset, length(cf))
    end

    # Add length assertion
    @assert length(gain) == length(cf)

    # Run model
    model!(
        stim,
        ffGn_hsr,
        ffGn_lsr,
        cf,
        n_chan,
        1 / fs,
        len_total,
        cohc,
        cihc,
        species_flag,
        100.0,
        powerlaw_mode,
        cn_tau_e,
        cn_tau_i,
        cn_delay,
        cn_amp,
        cn_inh,
        ic_tau_e,
        ic_tau_i,
        ic_delay,
        ic_amp,
        ic_inh,
        moc_cutoff,
        moc_beta,
        moc_offset,
        moc_minval,
        moc_maxval,
        moc_weight,
        moc_width,
        dur_pad_left,
        moc_delay,
        Int(moc_fix_gain),
        controlout,
        c1out,
        c2out,
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
        gainpostmix,
    )

    # Return
    outputs = [
        controlout,
        c1out,
        c2out,
        ihcout,
        expout_hsr,
        sout1_hsr,
        sout2_hsr,
        synout_hsr,
        expout_lsr,
        sout1_lsr,
        sout2_lsr,
        synout_lsr,
        hsrout,
        lsrout,
        cnout,
        icout,
        mocwdr,
        mocic,
        gain,
        gainpostmix
    ]
    if clip_left | clip_right
        outputs = map(outputs) do output
            output = map(output) do channel
                idx_left = clip_left ? (len_pad_left + 1) : 1
                idx_right = clip_right ? length(channel) - len_pad_right : length(channel)
                channel = channel[idx_left:idx_right]
            end
        end
    end
    return outputs
end

function sim_gfc2023!(mem::GFC2023_Mem, args...; fractional=false, clean=false, kwargs...)
    # If fractional, we should first fill the ffGn values in
    if fractional
        update_ffGn!(mem)
    end

    # If clean, we should zero out memory state before running
    if clean
        zero_state!(mem)
    end

    # Pass to sim_gfc2023! standard method
    return sim_gfc2023!(
        mem.ffGn_hsr,
        mem.ffGn_lsr,
        mem.controlout,
        mem.c1out,
        mem.c2out,
        mem.ihcout,
        mem.expout_hsr,
        mem.sout1_hsr,
        mem.sout2_hsr,
        mem.synout_hsr,
        mem.expout_lsr,
        mem.sout1_lsr,
        mem.sout2_lsr,
        mem.synout_lsr,
        mem.hsrout,
        mem.lsrout,
        mem.cnout,
        mem.icout,
        mem.mocwdr,
        mem.mocic,
        mem.gain,
        mem.gainpostmix,
        args...;
        fs=mem.fs,
        kwargs...
    )
end

"""
    sim_gfc2023_dict(args...; kwargs...)

Wrapper around `sim_gfc2023` that returns outputs as a dictionary.
"""
function sim_gfc2023_dict(args...; kwargs...)
    control, c1, c2, ihc, expon_hsr, sout1_hsr, sout2_hsr, syn_hsr, expon_lsr, sout1_lsr, sout2_lsr, syn_lsr, hsr, lsr, cn, ic, mocwdr, mocic, gain, gainpostmix = sim_gfc2023(args...; kwargs...)
    return Dict(
        "control" => control,
        "c1" => c1,
        "c2" => c2,
        "ihc" => ihc,
        "expon_hsr" => expon_hsr,
        "sout1_hsr" => sout1_hsr,
        "sout2_hsr" => sout2_hsr,
        "syn_hsr" => syn_hsr,
        "expon_lsr" => expon_lsr,
        "sout1_lsr" => sout1_lsr,
        "sout2_lsr" => sout2_lsr,
        "syn_lsr" => syn_lsr,
        "hsr" => hsr,
        "lsr" => lsr,
        "cn" => cn,
        "ic" => ic,
        "mocwdr" => mocwdr,
        "mocic" => mocic,
        "gain" => gain,
        "gainpostmix" => gainpostmix
    )
end

