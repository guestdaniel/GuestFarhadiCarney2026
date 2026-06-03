# GuestFarhadiCarney2026.jl
# Last updated: 06/02/2026
# Daniel R. Guest
#
# This package uses produces the figures and analyses reported in: Guest, D. R.,
# Farhadi, A., and Carney, L. H. (2026). "A computational model of the mammalian
# auditory periphery with a closed-loop medial olivocochlear reflex simulating
# across-channel efferent gain control." Journal of the Acoustical Society of
# America. In press.
#
# A NOTE ON MEMOIZATION:
# This script relies on an automated memoization system, see
# https://en.wikipedia.org/wiki/Memoization for more information about the
# concept. This can be seen in the @memo config macros before most simulation
# calls in the source code. The result is that the first run of this script may
# be very slow, as all computations have to be performed from scatch, but on
# subsequent runs everything will be virtually instantaneous (except for
# plotting time), as results can be automagically recalled from cache and need
# not be re-computed. 
# 
# A NOTE ON FIGURE SIZES:
# Plot sizes are initially specified in pixels, but implicitly we assume that
# because we will render at 300 DPI, that plot sizes are specified in
# inches*100... i.e., if we take plot size (600, 400)px and do px_per_unit=3 we
# achieve (1800, 1200)px and or (6, 4)in at 300 DPI. This is mostly handled by
# the `saveplot` function below.
#
# A NOTE ON JULIA:
# Even if you do not know Julia, this script should still be fairly
# approachable. Each individual figure is implemented in a function in src.jl or
# one of its child source files. If there are any functions that you wish to
# query, start by checking whether they have documentation by typing a ? mark in
# the REPL to enter help mode and then typing the name of the function. The
# @which and @less macros are also helpful. For a specific function call, you
# can use them to identify which method is used and the code that defines it,
#   e.g. `@which genfig_results_cas_physio_Warren1989a_Fig6()` 
#        `@less genfig_results_cas_physio_Warren1989a_Fig6()`

# Import packages 
using GuestFarhadiCarney2026
using Utilities
using DrWatson
using CairoMakie
using AuditorySignalUtils
using ColorSchemes

# Define project name and the `Config` settings (these determine behavior of 
# memoization features. The default behavior is to both save and load from 
# the file cache).
config = Config(; load_from_cache=true)

# Set the theme used for the figures 
set_theme!(; 
    fontsize=10.0,
    Axis=(
        rightspinevisible=false,
        topspinevisible=false,
        xgridvisible=false,
        ygridvisible=false,
    )
)

# Configure parallel pool
# Optionally, you can run the line below to run everything in parallel. The
# number passed to the @parallel macro controls how many additional processes we
# try to add. There is no value in adding more processes than there are threads
# available on your system. Also, each process requires some memory (usually
# ~300–400 MB), so available memory on your system may limit the number of
# processes you can add. 8 is a reasonable starting point for most modern
# workstations.
@parallel(8)

# ##############################################################################
# PRE-FIGURE: Choose parameters for the model
# Define sub-params
β = 0.045  # baseline β parameter, before CF adjustments 
moc_width = 0.9

# Create params structs used by subsequent functions
params = MOCParams4(; 
    moc_weight=x -> 1.0, 
    moc_beta=x -> β * peaknorm_gaussian(log2(x/3e3), 0.0, 2.5),
    moc_offset=x -> max(5.0 * log2(x/2e3) + 3.0, 3.0),  # slope previously 3
    moc_width=moc_width,
)

params_singlechannel = MOCParams4(; 
    moc_weight=x -> 1.0, 
    moc_beta=x -> β * peaknorm_gaussian(log2(x/3e3), 0.0, 2.5),
    moc_offset=x -> max(5.0 * log2(x/2e3) + 3.0, 3.0),  # slope previously 3
    moc_width=0.0,
)

params_flat = MOCParams4(; 
    moc_weight=x -> 1.0, 
    moc_beta=x -> β,
    moc_offset=x -> 2.5,
    moc_width=moc_width,
)

# ##############################################################################
# FIGURE 1: MODEL ARCHITECTURE AND PHILOSOPHY
# This figure shows the architecture and design of the model in a schematic
# format. This figure is generated in Inkscape and does not involve this script.
#
# Figure 1. (A) A schematic depiction of the conceptual model underlying the
# present simulation of the contralateral MOCR. On the contralateral side
# (left), an elicitor drives activity in AN fibers, which subsequently excite
# tonotopically similar cochlear nucleus neurons located on the same side as the
# elicitor and, after a midline cross, tonotopically similar MOC neurons on the
# opposite side as the elicitor.  These MOC neurons each then project to
# innervate multiple ipsilateral OHCs. Crucially, because of the sometimes wide
# tonotopic span of innervation patterns by individual MOC neurons (Brown,
# 2014), cochlear gain at an ipsilateral (right) recording site could be
# influenced by the activity of multiple MOC neurons tuned to different
# frequencies. (B) A series of schematics indicating contributions to MOC
# activity that could influence cochlear gain in the ipsilateral cochlea. From
# left to right: the ipsilateral MOCR, the contralateral MOCR, and descending
# inputs to MOC. Here, “electrode” indicates where an electrode would be placed
# to record from ipsilateral AN fibers, relative to the contralateral elicitor,
# as in Warren and Liberman (1989a). (C) Schematic depicting a single channel in
# the model with MOCR (c.f. Fig 1, Zilany and Bruce, 2006). Blue arrows indicate
# the path of the control signal for the MOCR (γ), which is used internally but
# also provided as an output for subsequent analysis or simulations. (D)
# Schematic depicting how CAS was approximated via two separate single-ear
# simulations (Simulation 1 and Simulation 2). First, a population response to
# the contralateral elicitor stimulus is simulated (Simulation 1, top). Then, a
# single-channel response to the ipsilateral probe stimulus is simulated at the
# CF of interest (Simulation 2, bottom). However, instead of allowing cochlear
# gain to be dynamically determined based on the MOCR (as in the single channel
# depicted in C), the MOCR control parameter (γ) is extracted from the outputs
# of Simulation 1 and used to determine cochlear gain at all time points during
# Simulation 2. Note that because this model architecture does not simulate any
# ipsilateral-to-contralateral pathways, the architecture omits any possible
# influence of ipsilateral responses on contralateral responses via the MOCR.

################################################################################
# FIGURE 2: RLF METHODS
# (A) Example RLF for a model HSR neuron. Gray arrows indicate rate threshold
# (TH) and RLF midpoint (MP), while gray brackets indicate the 10–90% rate range
# and the dynamic range, as defined in Sec II. (B) Example comparison of two
# RLFs simulated in a quiet context (black) or with a contralateral elicitor
# reducing cochlear gain (red). ΔL and ΔR are marked with horizontal and
# vertical gray arrows, respectively.

# TODO: migrate appropriate code here

# # ############################################################################
# FIGURE 3: MODEL INTERNALS
# (A) Stimulus waveforms for the CAS paradigm, with a probe tone presented at 5
# dB SPL (top) or an elicitor tone presented at 70 dB SPL (bottom). Note that
# vertical axes are not to scale. The vertical gray line here, and in B and C,
# indicates the onset of the elicitor. (B) Simulated ipsilateral AN response for
# an 8-kHz HSR fiber in response to the stimuli in A, either with the MOC
# pathway disabled (dashed gray) or the MOC pathway enabled (red). The fiber had
# a spontaneous rate just over 100 sp/s, which is visible in the response before
# probe onset. (C) Responses at intermediate model stages underlying MOC gain
# control in response to the stimuli in A. The top panel shows lowpass-filtered
# LSR rates from the 8-kHz contralateral CF, used as a proxy for sound-level
# information in the ascending auditory pathway. The middle panel shows the
# resulting “gain factor” \gamma in the same channel based on the MOC IO
# nonlinearity and following the spatial-smoothing step. The bottom panel shows
# the effective time-varying change in cochlear gain corresponding to each
# \gamma value in time, as inferred from midpoint shifts in on-CF RLFs (see Sec.
# IIC for details). (D) The IO nonlinearity relating MOC rate to gain factor at
# the 8-kHz CF. Note that this subfigure shows the single-channel gain factor
# before the spatial-smoothing step; in practice, the smoothing step
# substantially attenuated final gain factors, at least for pure-tone
# stimulation for which stimulation was highly tonotopically specific. (E) The
# relationship between gain factor and effective change in cochlear gain (ΔL),
# as mentioned in C. (F) Parameters in the MOC model as a function of CF, either
# for the version with the same parameters at all CFs (gray dashed) or
# CF-varying parameters (black). The IO threshold parameter (labeled “Offset”)
# is shown on the top and the IO slope parameter (labelled “\beta”) on the
# bottom.

# >>>>> FIGURE 2A/B/C: Response examples
genfig_results_intro_example_responses(params; config=config)

# >>>>> FIGURE 2D/E/F: Input-output relationships
genfig_results_intro_io_stack(params, params_flat; config=config)


# ##############################################################################
# FIGURE 4: CAS THRESHOLDS VS WARREN AND LIBERMAN (1989)
# (A) Example iso-response tuning curve for a 2-kHz model HSR fiber (black line)
# in comparison to the frequency and level of a probe stimulus (gray marker) and
# frequency and level range of an elicitor stimulus (red line) used to determine
# the CAS threshold. (B) Corresponding firing rate versus elicitor level curve
# for the model neuron in A, with a red arrow indicating CAS threshold (elicitor
# level achieving 5% reduction in probe rate). (C) Effect of the MOC IO
# nonlinearity parameter \beta, relative to the default value, on the curve
# shown in B. (D) CAS thresholds as function of CF from Warren and Liberman
# (1989a), Fig. 6. Colored markers indicate individual-fiber data; the solid
# gray line a trendline minimizing mean-squared-error across all fiber types;
# the dashed gray line a “lower envelope” of the data. (E) Data from A replotted
# in gray alongside model LSR thresholds as a function of CF (dashed purple
# line) and CAS thresholds for a model HSR fiber (solid red line).

# >>>>> FIGURE 3A: Example tuning curve and stimulus schematic
genfig_results_cas_sim_tc(; config=config, savefig=false)

# >>>>> FIGURE 3B: Example SLF
genfig_results_cas_sim_slf(params; config=config, savefig=false)

# >>>>> FIGURE 3C: SLF vs params
genfig_results_cas_sim_slf_vs_params(params; config=config, savefig=false)

# >>>>> FIGURE 3D: Physiological data
genfig_results_cas_physio_Warren1989a_Fig6(; savefig=false)

# >>>>> FIGURE 3E: Simulated data
genfig_results_cas_threshold_vs_cf(
    [params_flat, params]; 
    config=config, 
    savefig=false,
)

# ##############################################################################
# FIGURE 5: CAS THRESHOLDS VS WARREN AND LIBERMAN (1989)
# (A) Example tuning curve for a 2-kHz model HSR fiber (black line) in
# comparison to the frequency and level of a probe stimulus (gray line) and
# frequency and level range of an elicitor stimulus (red marker) used to
# determine ΔL and ΔR from an RLF simulation. (B) ΔL for elicitors in the range
# of 60–85 dB SPL from Figure 2 of Warren and Liberman (1989b). Markers indicate
# individual-fiber ΔLs. The solid gray line is a trendline minimizing
# mean-squared-error across all fiber types. The dashed gray line is an “upper
# envelope” of the data. Markers in red show simulated ΔLs. (C) ΔR for elicitors
# in the range of 60–70 dB SPL from Figure 7 of Warren and Liberman (1989a).
# Markers indicate individual-fiber ΔRs. The solid gray curve is a trendline
# minimizing mean-squared-error across all fiber types. The dashed gray curve is
# an “upper envelope” of the data. The horizontal gray line simply denotes 0%
# ΔR. Markers in red show simulated ΔRs.

# >>>>> FIGURE 4A: Example tuning curve and stimulus schematic
genfig_results_cas_sim_tc(; config=config, fn="fig_results_cas_sim_tc_2.png", annotate_flag=2, size=(2.5*100, 2*100), xlims=(0.1, 30.0))

# >>>>> FIGURE 4B: RATE-LEVEL FUNCTION
genfig_results_rlf(; config=config)

# >>>>> FIGURE 4C: RLF SHIFT due to CAS
genfig_results_intro_rlf_shift(params; config=config)

# >>>>> FIGURE 4D: ΔL vs CF 
genfig_results_cas_ΔL_vs_cf([params_flat, params]; config=config, savefig=true)

# >>>>> FIGURE 4E: MAGNITUDE VS CF
genfig_results_cas_magnitude_vs_cf([params_flat, params]; config=config, savefig=true)

# ##############################################################################
# FIGURE 6: ALL WIDEBAND-RELATED FIGURES
# (A) Example tuning curve for a 2-kHz model HSR fiber (black line) in
# comparison to the frequency and level of a probe stimulus (gray marker) and
# level and frequency range of an elicitor stimulus (red line and markers) used
# to characterize the tuning characteristics of MOC gain control in the model.
# (B) Reproduction of an example LSR fiber dataset from Warren and Liberman
# (1989b), Figure 6. (C) Simulation of results from B using a 2-kHz LSR fiber.
# Probe level matched the midpoint of the 2-kHz RLF while elicitor level was
# fixed at 70 dB SPL. (D) Effect of the MOC smoothing bandwidth, expressed in
# units of standard deviation in octaves, on result shown in C but with an HSR
# fiber, with different values indicated by different colors. Here, the y-axis
# was started at 75 sp/s rather than 0 sp/s to aid comparison to the LSR fibers
# shown in B and C, where the bottom of the axis (0 sp/s) and spontaneous rate
# (0.6 sp/s for B or ~1 sp/s for C) are similar. (E) Scatter plot of the best
# suppressor frequency (BSF) vs CF for a range of CFs, given a 65 dB SPL
# pure-tone suppressor. The straight line indicates where points would fall if
# BSF were always equal to CF. (F) Effect of widening the bandwidth of a
# bandlimited noise suppressor on the result shown in C. Colors are as in D.
# Lines show a loess interpolation of the firing rates in response to noise,
# separately for each tested parameter value.

# Select parameter sets
params_varying_width = map([0.0, 0.2, 0.4, 0.8, 1.6]) do width
    MOCParams4(
        moc_weight=x -> 1.0, 
        moc_beta=x -> β * peaknorm_gaussian(log2(x/3e3), 0.0, 2.5),
        moc_offset=x -> max(5.0 * log2(x/2e3) + 3.0, 3.0),  # slope previously 3
        moc_width=width,
    )
end

# >>>>> FIGURE 5A: Example tuning curve and stimulus schematic
genfig_results_cas_sim_tc(; 
    config=config, 
    fn="fig_results_cas_sim_tc_3.png",
    annotate_flag=3,
)

# >>>>> FIGURE 5B: Physiological data example
fig, ax = plot_Warren1989b_Fig6(;
    size=(2.0*100, 1.6*100),
)
saveplot("fig_results_cas_physio_Warren1989b_Fig6.png", fig)

# >>>>> Figure 5C: Physiological data comparison
genfig_results_cas_sim_freq_sweeps(params; config=config, savefig=true)

# >>>>> Figure 5D: Effect of parameters 
genfig_results_cas_sim_freq_sweeps(params_varying_width; config=config, savefig=true)

# >>>>> Figure 5E: CF vs BSF
genfig_results_cas_sim_cf_vs_bsf([params]; config=config, savefig=true)

# >>>>> Figure 5F: Band widening experiment
genfig_results_cas_band_widening(params_varying_width; config=config, savefig=true)

# ##############################################################################
# SUPP FIG 1: LEVEL GROWTH AND RATE LEVEL FUNCTIONS
# CAPTION
# ##############################################################################
genfig_supp_level_growth_panel(; config=config)

genfig_supp_level_growth_panel(; 
    config=config, 
    cfs=LogRange(3.0e3, 6.0e3, 41), 
    xlims=(3.0, 6.0), 
    highlight_inset=false, 
    xticks=2.0:1.0:6.0, 
    xminorticksvisible=false,
    fig=Figure(; size=(650, 200)),
    ylabels=fill("", 4),
    fn="fig_suppelmental_growth_curves_inset.png"
)

# ##############################################################################
# FIGURE 6: Pretty super figure summarizing overall MOC behavior across elicitor level
# *Sketch*
params = fetch_params("best")

# Rate-level function rainbow
genfig_results_summary_rlf(params; config=config)

# Plot maxium attenuations
genfig

# ##############################################################################
# FOLLOW-UP FIGURE: Measure average rate responses for elicitor tones in contralateral ear 
# Fetch model
genfig_maximum_attenuation(params; config=config)