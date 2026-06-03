# This script tests whether the power-law adaptation approximation code described in Guest
# and Carney (2024) can truly approximate the desired power-law adaptation. We compare 2014
# model responses under true power-law adaptation in several stages to new model responses
# using the approximate scheme with optimized weights. A few stages are not exactly matched
# and are excluded from the test. The test stimulus is a mid-level pure tone. 
@testset "Exponential process approximation to PLA" begin
    @testset "Sampling rate = $fs Hz" for fs in [100e3]
        cfs = [500.0, 1000.0, 2000.0, 4000.0]
        @testset "CF: $cf" for cf in cfs
            orig, new = run_2014_vs_2023(
                pt(cf, 50.0, 0.2, fs), 
                [cf],
                Dict{Symbol, Any}(),
                Dict{Symbol, Any}(:powerlaw_mode => 2, :fs => fs, :moc_weight => 0.0, :dur_pad_left => 0.2),
            )
            @testset "$cf-kHz pure tone, 50 dB SPL, stage: $stage" for stage in testparams["stages_peripheral"] 
                @test isapprox(orig[stage][1], new[stage][1]; rtol=0.03) broken=( ((cf==4e3) & ((stage == "lsr") | (stage == "hsr"))) | ((cf==0.5e3) & ((stage == "c1") | (stage == "c2"))) )
            end
        end
    end
end
