@testset "Single- vs multi-channel" begin
    cfs = [500.0, 1000.0, 2000.0, 4000.0]
    stim = pt(1000.0, 50.0)
    multichannel = sim_gfc2023_dict(pt(1000.0, 50.0), cfs; moc_weight=0.0)
    @testset "CF = $cf Hz" for (idx_cf, cf) in enumerate(cfs)
        singlechannel = sim_gfc2023_dict(stim, [cf]; moc_weight=0.0)
        @testset "stage: $stage" for stage in testparams["stages_peripheral"] 
            @test isapprox(multichannel[stage][idx_cf], singlechannel[stage][1]; rtol=0.001)
        end
    end
end
