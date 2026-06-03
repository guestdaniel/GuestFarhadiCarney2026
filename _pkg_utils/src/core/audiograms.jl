export fit_audiogram, Audiogram, string, rms, NH, N1, N2, N3, N4, N5, N6, N7, S1, S2, S3, F1, F2, F3, F4, interpolate

if Sys.iswindows()
    data_cat = matread("C:\\Users\\dguest2\\cl_code\\urear_2020b\\thresholds_interpolated_cat.mat")
    data_human_shera = matread("C:\\Users\\dguest2\\cl_code\\urear_2020b\\thresholds_interpolated_human_shera.mat")
    data_human_glasberg_moore = matread("C:\\Users\\dguest2\\cl_code\\urear_2020b\\thresholds_interpolated_human_glasberg_moore.mat")
else
    data_cat = matread("/home/dguest2/thresholds_interpolated_cat.mat")
    data_human_shera = matread("/home/dguest2/thresholds_interpolated_human_shera.mat")
    data_human_glasberg_moore = matread("/home/dguest2/thresholds_interpolated_human_glasberg_moore.mat")
end

"""
    fit_audiogram(freqs, losses, species[, ohc_loss_target])

# Arguments
- `freqs`: Array of audiogram frequencies (Hz)
- `losses`: Array of audiometric threshold shifts (dB re: 0 dB HL)
- `species`: String indicating species of the AN model ("cat", "human", or "human_glasberg")
- `ohc_loss_target`: Threshold shift attributable to OHC loss alone
"""
function fit_audiogram(
    freqs::AbstractVector{Float64},
    losses::AbstractVector{Float64},
    species::String,
    ohc_loss_target::AbstractVector{Float64}=(2/3) .* losses,
)
    # Load data relating IHC and OHC loss to thresholds
    if species == "cat"
        data = data_cat
    elseif species == "human"
        data = data_human_shera
    else
        data = data_human_glasberg_moore
    end
    cf = dropdims(data["cf_out"]; dims=1)
    cihc = dropdims(data["cihc_out"]; dims=1)
    cohc = dropdims(data["cohc_out"]; dims=1)
    threshold = data["threshold_out"]

    # Express thresholds in terms of shift re: absolute threshold
    threshold_shift = threshold .- threshold[:, 1, 1]

    # Pre-allocate storage
    cohc_out = zeros(length(freqs))
    cihc_out = zeros(length(freqs))
    loss_ihc = zeros(length(freqs))
    loss_ohc = zeros(length(freqs))

    for (idx_freq, freq) in enumerate(freqs)
        # Identify index into CF nearest to current frequency
        idx_cf = argmin(abs.(cf .- freq))

        # Handle edge case where requested OHC-induced loss is greater than possible
        if ohc_loss_target[idx_freq] > threshold_shift[idx_cf, 1, end]
            cohc_out[idx_freq] = 0.0
            loss_ohc[idx_freq] = threshold_shift[idx_cf, 1, end]
        # Otherwise, set OHC loss based on requested OHC-induced threshold shift
        else
            idx = argmin(abs.(threshold_shift[idx_cf, 1, :] .- ohc_loss_target[idx_freq]))
            cohc_out[idx_freq] = cohc[idx]
            loss_ohc[idx_freq] = threshold_shift[idx_cf, 1, idx]
        end

        # Set IHC-induced loss to account for remaining threshold shift above OHC-induced threshold shift
        idx = argmin(abs.(cohc .- cohc_out[idx_freq]))
        loss_ihc[idx_freq] = losses[idx_freq] - loss_ohc[idx_freq]

        # Handle edge case (as above)
        if losses[idx_freq] > threshold_shift[idx_cf, end, idx]
            cihc_out[idx_freq] = 0.0
        else
        # Set IHC loss
            idx2 = argmin(abs.(threshold_shift[idx_cf, :, idx] .- losses[idx_freq]))
            cihc_out[idx_freq] = cihc[idx2]
        end
    end

    return cohc_out, cihc_out
end

# Audiogram type
@with_kw struct Audiogram <: Component
    freqs::Vector{Float64}=[250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0]
    thresholds::Vector{Float64}=zeros(length(freqs))
    species::String="human"
    desc::String="NH"
end

Base.string(audiogram::Audiogram) = audiogram.desc

function Audiogram(thresholds::Vector{Float64}, species::String="human", desc::String="none")
    Audiogram([250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0], thresholds, species, desc)
end

function fit_audiogram(audiogram::Audiogram, cfs)
    itp = LinearInterpolation(audiogram.freqs, audiogram.thresholds; extrapolation_bc=Interpolations.Flat())
    losses = itp.(cfs)
    fit_audiogram(cfs, losses, audiogram.species)
end

# Canonical audiograms
aud_freqs = [250.0, 375.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 6000.0]
NH = Audiogram(aud_freqs, zeros(length(aud_freqs)), "human", "NH")
N1 = Audiogram(aud_freqs, [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 15.0, 20.0, 30.0, 40.0], "human", "HI-N1")
N2 = Audiogram(aud_freqs, [20.0, 20.0, 20.0, 22.5, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0], "human", "HI-N2")
N3 = Audiogram(aud_freqs, [35.0, 35.0, 35.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0], "human", "HI-N3")
N4 = Audiogram(aud_freqs, [55.0, 55.0, 55.0, 55.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0], "human", "HI-N4")
N5 = Audiogram(aud_freqs, [65.0, 67.5, 70.0, 72.5, 75.0, 80.0, 80.0, 80.0, 80.0, 80.0], "human", "HI-N5")
N6 = Audiogram(aud_freqs, [75.0, 77.5, 80.0, 82.5, 85.0, 90.0, 90.0, 95.0, 100.0, 100.0], "human", "HI-N6")
N7 = Audiogram(aud_freqs, [90.0, 92.5, 95.0, 100.0, 105.0, 105.0, 105.0, 105.0, 105.0, 105.0], "human", "HI-N7")
S1 = Audiogram(aud_freqs, [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 15.0, 30.0, 55.0, 70.0], "human", "HI-S1")
S2 = Audiogram(aud_freqs, [20.0, 20.0, 20.0, 22.5, 25.0, 35.0, 55.0, 75.0, 95.0, 95.0], "human", "HI-S2")
S3 = Audiogram(aud_freqs, [30.0, 30.0, 35.0, 47.5, 60.0, 70.0, 75.0, 80.0, 80.0, 85.0], "human", "HI-S3")
F1 = Audiogram(aud_freqs, 20.0 .* ones(length(aud_freqs)), "human", "HI-F1")
F2 = Audiogram(aud_freqs, 40.0 .* ones(length(aud_freqs)), "human", "HI-F2")
F3 = Audiogram(aud_freqs, 60.0 .* ones(length(aud_freqs)), "human", "HI-F3")
F4 = Audiogram(aud_freqs, 80.0 .* ones(length(aud_freqs)), "human", "HI-F4")

# Functions to linearly interpolate between two audiograms, or series of audiograms
function interpolate(a1::Audiogram, a2::Audiogram)
    # First, verify that the two audiograms have the same underlying sampling
    @assert all(a1.freqs .== a2.freqs)

    # Next, take an average of the threshold at each frequency
    Audiogram(
        a1.freqs,
        map(mean, zip(a1.thresholds, a2.thresholds)),
        a1.species,
        a1.desc * " ➡  " * a2.desc,
    )
end

function interpolate(a1::Audiogram, a2::Audiogram, n_step::Int)
    # First, verify that the two audiograms have the same underlying sampling
    @assert all(a1.freqs .== a2.freqs)

    # Construct series of audiogram
    interpolated_audiograms = map(1:n_step) do idx_step
        # Get frequencies and species
        freqs = a1.freqs
        species = "human"

        # Linearly interpolate based on which step we're on
        hl = map(enumerate(freqs)) do (idx_freq, freq)
            a1.thresholds[idx_freq] + idx_step * (a2.thresholds[idx_freq]-a1.thresholds[idx_freq])/(n_step+1)
        end

        # Get description using SHA256 and return
        desc = bytes2hex(sha256(string(freqs) * string(hl) * string(species)))
        Audiogram(freqs, hl, species, desc)
    end

    # Return
    vcat(a1, interpolated_audiograms..., a2)
end

function interpolate(audiograms::Vector{Audiogram}, n_step::Int)
    temp = map(1:length(audiograms)-1) do idx
        if idx == 1
            interpolate(audiograms[idx], audiograms[idx+1], n_step)
        else
            interpolate(audiograms[idx], audiograms[idx+1], n_step)[2:end]
        end
    end
    vcat(temp...)
end