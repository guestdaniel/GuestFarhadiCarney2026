using Test
using AuditorySignalUtils
using ZilanyBruceCarney2014
using AuditoryMidbrain
using Statistics
using DSP
using Helios

# Test implementation of SFIE in C
#include("sfie_implementation.jl")

# Test regression against 2014 model in single-channel simulations
include("regression_2014_single_channel.jl")

# Test that multichannel model results are identical to single-channel model results
include("singlechannel_vs_multichannel.jl")

# Test differences between true and approximate power-law adaptation
include("approximate_powerlaw.jl")

# Test implementation of normal PDF function in C
include("normal_pdf.jl")

# ==========================================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Check whether new model outputs match 2014 model outputs (multichannel)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ==========================================================================================
# @testset "Regression vs 2014 --- multichannel" begin
#     # ======================================================================================
#     # Check response to 1 kHz pure tone at 1 kHz and 2 kHz CFs
#     # ======================================================================================
#     @testset "1-kHz pure tone, 50 dB SPL, multichannel, stage: $stage" for stage in stages_peripheral
#         orig, new = run_2014_vs_2023_pure_tone([1000.0, 2000.0], stage)
#         pairs = zip(orig, new)
#         @test all(map(pair -> isapprox(pair[1], pair[2]; rtol=get_rtol(stage)), pairs))
#     end
# end

# # ==========================================================================================
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~ Check whether subcortial model outputs look reasonable and are behaving
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ==========================================================================================
# @testset "Regression vs 2004 --- single channel" begin
#     # ======================================================================================
#     # Check response to 1 kHz pure tone at cochlear nucleus
#     # ======================================================================================
#     @testset "CN at parameter values: $params" for params in params_sfie
#         # Extract params
#         τ_e, τ_i, d, a, s = params

#         # Create stimulus
#         x = pt(1000.0, 50.0)

#         # Simulate response from C subcortical model
#         out = sim_gfc2023_dict(
#             x, 
#             1000.0;
#             cn_tau_e=τ_e,
#             cn_tau_i=τ_i,
#             cn_delay=d,
#             cn_amp=a,
#             cn_inh=s,
#             dur_pad_left=0.0,
#             clip_left=false,
#         )
#         new = out["cn"]

#         # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
#         old = sim_sfie_nc2004(
#             out["hsr"], 
#             τ_e=τ_e,
#             τ_i=τ_i,
#             d_i=d,
#             S=s,
#             A=a,
#         )

#         # Compare
#         @test old ≈ new
#     end
#     # ======================================================================================
#     # Check response to 1 kHz pure tone at inferior colliculus
#     # ======================================================================================
#     @testset "IC at parameter values: $params" for params in params_sfie
#         # Extract params
#         τ_e, τ_i, d, a, s = params

#         # Create stimulus
#         x = pt(1000.0, 50.0)

#         # Simulate response from C subcortical model
#         out = sim_gfc2023_dict(
#             x, 
#             1000.0;
#             ic_tau_e=τ_e,
#             ic_tau_i=τ_i,
#             ic_delay=d,
#             ic_amp=a,
#             ic_inh=s,
#             dur_pad_left=0.0,
#             clip_left=false,
#         )
#         new = out["ic"]

#         # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
#         old = sim_sfie_nc2004(
#             out["hsr"], 
#             τ_e=0.5e-3,
#             τ_i=2.0e-3,
#             d_i=1.0e-3,
#             S=0.6,
#             A=1.5,
#         )
#         old = sim_sfie_nc2004(
#             old,
#             τ_e=τ_e,
#             τ_i=τ_i,
#             d_i=d,
#             S=s,
#             A=a,
#         )

#         # Compare
#         @test old ≈ new
#     end
# end

# @testset "Regression vs 2004 --- multichannel" begin
#     # ======================================================================================
#     # Check response to 1 kHz pure tone at inferior colliculus
#     # ======================================================================================
#     @testset "CN at parameter values: $params" for params in params_sfie
#         # Extract params
#         τ_e, τ_i, d, a, s = params

#         # Create stimulus
#         x = pt(1000.0, 50.0)

#         # Simulate response from C subcortical model
#         out = sim_gfc2023_dict(
#             x, 
#             [1000.0, 2000.0];
#             ic_tau_e=τ_e,
#             ic_tau_i=τ_i,
#             ic_delay=d,
#             ic_amp=a,
#             ic_inh=s,
#             dur_pad_left=0.0,
#             clip_left=false,
#         )
#         new = out["ic"]

#         # Simulate response from AuditoryMidbrain.jl for cochlear nucleus stage
#         old = map(out["hsr"]) do x 
#             cn = sim_sfie_nc2004(
#                 x, 
#                 τ_e=0.5e-3,
#                 τ_i=2.0e-3,
#                 d_i=1.0e-3,
#                 S=0.6,
#                 A=1.5,
#             )
#             ic = sim_sfie_nc2004(
#                 cn,
#                 τ_e=τ_e,
#                 τ_i=τ_i,
#                 d_i=d,
#                 S=s,
#                 A=a,
#             )
#             return ic
#         end

#         # Compare old to new
#         pairs = zip(old, new)
#         @test all(map(pair -> isapprox(pair[1], pair[2]), pairs))
#     end
# end

# # ==========================================================================================
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~ Check whether Julia and Mex wrappers provide same outputs
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ==========================================================================================
# # @testset "Julia vs Mex" begin
# #     # First, we'll compare the Julia wrapper to the Mex wrapper by simulating responses to a
# #     # pure tone with gain control disabled and verifying that responses at each output stage
# #     # available in the Mex wrapper (ihc, hsr, lsr, ic, and gain) produce matched outputs for
# #     # a single-CF response
# #     @testset "1-kHz pure tone, 50 dB SPL, gain control disabled, single channel, stage: $stage" for stage in ["ihc", "hsr", "lsr", "ic", "gain"]
# #         # Synthesize pure tone
# #         x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.3, 100e3), 0.01, 100e3), 50.0)

# #         # Run both models with gain control disabled
# #         julia = sim_gfc2023_dict(
# #             x, 
# #             1000.0; 
# #             dur_pad_left=0.0, 
# #             dur_pad_right=0.0,
# #             moc_weight_ic=0.0,
# #             moc_weight_wdr=0.0,
# #         )[stage]
# #         matlab = sim_gfc2023_wrapper_dict_mex(
# #             x, 
# #             1000.0;
# #             moc_weight_ic=0.0,
# #             moc_weight_wdr=0.0,
# #         )[stage]
# #         @test isapprox(julia, matlab; rtol=get_rtol(stage))
# #     end

# #     # Next, we'll repeat the same simulations above except that we will turn gain control 
# #     # on with very typical parameter values
# #     @testset "1-kHz pure tone, 50 dB SPL, gain control enabled, single channel, stage: $stage" for stage in ["ihc", "hsr", "lsr", "ic", "gain"]
# #         # Synthesize pure tone
# #         x = scale_dbspl(cosine_ramp(pure_tone(1000.0, 0.0, 0.3, 100e3), 0.01, 100e3), 50.0)

# #         # Run both models with gain control disabled
# #         julia = sim_gfc2023_dict(
# #             x, 
# #             1000.0; 
# #             dur_pad_left=0.0, 
# #             dur_pad_right=0.0,
# #             moc_weight_ic=1.0,
# #             moc_weight_wdr=1.0,
# #         )[stage]
# #         matlab = sim_gfc2023_wrapper_dict_mex(
# #             x, 
# #             1000.0;
# #             moc_weight_ic=1.0,
# #             moc_weight_wdr=1.0,
# #         )[stage]
# #         @test isapprox(julia, matlab; rtol=get_rtol(stage))
# #     end
# # end

