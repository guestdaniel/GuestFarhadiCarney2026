export sim_orig, sim_orig_dict

function sim_orig(
    x::Vector{Float64}, 
    cf::Float64; 
    fs::Float64=100e3,
    cohc::Float64=1.0,
    cihc::Float64=1.0,
    species::String="human",
    power_law::String="actual", 
    fractional::Bool=false,
)
    # Convert human-readable arguments into C-side floats/ints
    species_flag = Dict(
        "cat" => 1,
        "human" => 2,
        "human_glasberg" => 3
    )[species]

    implnt = Dict(
        "actual" => 1.0,
        "approximate" => 0.0
    )[power_law]

    noiseType = Dict(
        true => 1.0,
        false => 0.0
    )[fractional]

    # Synthesize ffGn
    if noiseType == 1.0
        ffGn = ffGn_native(
            Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)),
            1/fs_synapse,
            0.9,
            noiseType,
            100.0,
        )
    else
        ffGn = zeros(Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)))
    end

    # Call original model functions
    ihcout = zeros(length(x))
    c1out = zeros(length(x))
    c1vihcout = zeros(length(x))
    c2out = zeros(length(x))
    c2vihcout = zeros(length(x))
    controlout = zeros(length(x))
    ccall(
        (:IHCDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Cdouble,             # cf
            Cint,                # nrep
            Cdouble,             # tdres
            Cint,                # totalstim
            Cdouble,             # cohc
            Cdouble,             # cihc
            Cint,                # species
            Ptr{Cdouble},        # ihcout
            Ptr{Cdouble},        # c1out
            Ptr{Cdouble},        # c1vihcout
            Ptr{Cdouble},        # c1out
            Ptr{Cdouble},        # c2vihcout
            Ptr{Cdouble},        # controlout
        ),
        x, cf, 1, 1/fs, length(x), cohc, cihc, species_flag, ihcout, c1out, c1vihcout, c2out, c2vihcout, controlout, # pass arguments
    )

    synout = zeros(length(ihcout))
    exponout = zeros(length(x))
    delaypoint = Int(floor(7500 / (cf / 1e3)))
    powerlawin = zeros(length(x) + delaypoint*3)
    sout1 = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    sout2 = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    len_noise = Int(ceil((length(ihcout) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3))
    ffGn = zeros(len_noise)
    ccall(
        (:SYNAPSEDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Ptr{Cdouble},        # randNums
            Cdouble,             # tdres
            Cdouble,             # cf
            Cint,                # totalstim
            Cint,                # nrep
            Cdouble,             # spont
            Cdouble,             # noisetype
            Cdouble,             # implementation
            Cdouble,             # sampFreq
            Ptr{Cdouble},        # synout
            Ptr{Cdouble},        # exponout
            Ptr{Cdouble},        # powerlawin
            Ptr{Cdouble},        # sout1
            Ptr{Cdouble},        # sout2
            Ptr{Cvoid},          # decimate function handle
        ),
        ihcout, ffGn, 1/fs, cf, length(ihcout), 1, 100.0, noiseType, implnt, 10e3, synout, exponout, powerlawin, sout1, sout2, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
    )

    if noiseType == 1.0
        ffGn = ffGn_native(
            Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)),
            1/fs_synapse,
            0.9,
            noiseType,
            0.1,
        )
    else
        ffGn = zeros(Int(ceil((length(x) + 2 * floor(7500 / (cf / 1e3))) * 1/fs * 10e3)))
    end
    synout_lsr = zeros(length(ihcout))
    exponout_lsr = zeros(length(x))
    powerlawin_lsr = zeros(length(x) + delaypoint*3)
    sout1_lsr = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    sout2_lsr = zeros(Int(ceil((length(ihcout)+2*delaypoint) * 1/100e3 * 10e3)))
    ccall(
        (:SYNAPSEDEBUG, "C:\\Users\\dguest2\\cl_code\\Helios\\external\\julia\\libzbc2014debug.so"),
        Cvoid,                   # return type
        (                        # arg types
            Ptr{Cdouble},        # px
            Ptr{Cdouble},        # randNums
            Cdouble,             # tdres
            Cdouble,             # cf
            Cint,                # totalstim
            Cint,                # nrep
            Cdouble,             # spont
            Cdouble,             # noisetype
            Cdouble,             # implementation
            Cdouble,             # sampFreq
            Ptr{Cdouble},        # synout
            Ptr{Cdouble},        # exponout
            Ptr{Cdouble},        # powerlawin
            Ptr{Cdouble},        # sout1
            Ptr{Cdouble},        # sout2
            Ptr{Cvoid},          # decimate function handle
        ),
        ihcout, ffGn, 1/fs, cf, length(ihcout), 1, 0.1, noiseType, implnt, 10e3, synout_lsr, exponout_lsr, powerlawin_lsr, sout1_lsr, sout2_lsr, @cfunction(decimate, Ptr{Cdouble}, (Ptr{Cdouble}, Cint, Cint)),
    )

    hsr = synout ./ (1.0 .+ 0.75e-3 .* synout)
    lsr = synout_lsr ./ (1.0 .+ 0.75e-3 .* synout_lsr)

    # Return
    return controlout, c1out, c1vihcout, c2out, c2vihcout, ihcout, synout, exponout, powerlawin, sout1, sout2, hsr, lsr
end

function sim_orig_dict(args...; kwargs...)
    control, c1, c1vihc, c2, c2vihc, ihc, syn, expon, powerlaw, sout1, sout2, hsr, lsr = sim_orig(args..., kwargs...)
    return Dict(
        "control" => control,
        "c1" => c1,
        "c1vihc" => c1vihc,
        "c2" => c2,
        "c2vihc" => c2vihc,
        "ihc" => ihc,
        "syn" => syn,
        "expon" => expon,
        "powerlaw" => powerlaw,
        "sout1" => sout1,
        "sout2" => sout2,
        "hsr" => hsr,
        "lsr" => lsr,
    )
end
