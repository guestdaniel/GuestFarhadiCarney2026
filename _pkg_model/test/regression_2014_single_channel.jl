# This script tests whether the current model code correctly matches responses from the old
# model when parameters are configured accordingly. There are several points of awkwardness
# that arise when comparing model responses from the different stages, these are mostly
# handled by the parameter settings or by the model-response postprocessing code. A few
# stages with known issue are excluded. True power-law adaptation is used by the new model;
# testing for the power-law approximation provided by `powerlaw_mode=2` is tested elsewhere. 
@testset "Regression vs 2014 --- single channel" begin
    # Loop over different frequencies
    @testset "CF/tone freq = $cf Hz" for cf in [1000.0, 2000.0, 4000.0]
        # Simulate responses
        orig, new = run_2014_vs_2023(
            pt(cf, 50.0), 
            [cf],
            Dict{Symbol, Any}(),  # leave old param values at default values
            Dict{Symbol, Any}(:powerlaw_mode => 1, :moc_weight => 0.0, :dur_pad_left => 0.2),
        )

        # Loop through each stage and verify match
        @testset "stage: $stage" for stage in testparams["stages_peripheral"] 
            # If the stage is one of [sout1, sout2, syn, expon], we need to handle the name
            # accordingly; for the efferent model, we have stage_hsr and stage_lsr, but for
            # the 2014 model we only have stage...
            if occursin("_", stage)
                stage_2014 = split(stage, "_")[1]
                stage_2023 = stage 
            else
                stage_2014 = stage
                stage_2023 = stage 
            end
            @test isapprox(orig[stage_2014][1], new[stage_2023][1]; rtol=rtol_2014_vs_2023(stage)) broken=((cf==4e3) & (stage=="lsr"))
        end
    end
end
