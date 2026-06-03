using Helios
using UnicodePlots
function ihcnl(x) 
    ccall(
        (:NLogarithm, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"),
        Cdouble,
        (
            Cdouble,
            Cdouble,
            Cdouble,
            Cdouble,
        ),
        x,
        0.1,
        3.0,
        1000.0,
    )
end

x = LinRange(-0.001, 0.001, 1000)
lineplot(x, ihcnl.(x))