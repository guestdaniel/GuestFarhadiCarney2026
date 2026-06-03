using Test
using Helios
using Statistics
using Distributions
using Libdl

@testset "Normal PDF" begin
    @testset "Params: $param" for param in [(0.0, 1.0), (-2.0, 2.0), (1.0, 0.01)]
        # Select input range
        x = -8.0:0.1:8.0

        # Calculate using Julia libraries
        y = pdf(Normal(param...), x)

        # Calculate using direct implementation in C
        lib = Libdl.dlopen("C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so")
        modelfunc = Libdl.dlsym(lib, :normal_pdf)
        ŷ = map(x) do _x
            ccall(
                modelfunc,
                Cdouble,
                (Cdouble, Cdouble, Cdouble),
                _x,
                param[1],
                param[2],
            )
        end
        Libdl.dlclose(lib)

        @test all(y .≈ ŷ)
    end
end
