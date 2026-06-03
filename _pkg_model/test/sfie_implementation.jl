@testset "SFIE implementation" begin
    # ======================================================================================
    # Check that get_alpha_norm produces identical values in C implementation
    # ======================================================================================
    @testset "Computing filter coefficients for τ = $τ" for τ in [1e-3, 2e-3, 4e-3]
        @test begin
            b_julia, a_julia = AuditoryMidbrain.get_α_normalized(τ, 100e3, 1.0)
            b_c = zeros(2); a_c = zeros(3);
            ccall(
                (:get_alpha_norm, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Cdouble,
                    Cdouble,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                τ, 100e3, 1.0, b_c, a_c,
            )
            (b_c ≈ b_julia) & (a_c ≈ a_julia)
        end
    end

    # ======================================================================================
    # Check that filter_alpha produces identical values in C implementation
    # ======================================================================================
    @test begin
        # Julia version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_julia = filt(b, a, x)

        # C version
        b, a = AuditoryMidbrain.get_α_normalized(1e-3, 100e3, 1.0)
        x = sin.(2π .* 250.0 .* (0.0:(1/100e3):(0.1 - 1/100e3)))
        y_c = zeros(length(x))
        for i = 0:1:(length(x)-1)
            ccall(
                (:filter_alpha, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
                Cvoid,
                (
                    Ptr{Cdouble},
                    Cint,
                    Cdouble,
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                    Ptr{Cdouble},
                ),
                x, i, 100e3, b, a, y_c,
            )
        end

        # Compare
        y_c ≈ y_julia
    end
end
