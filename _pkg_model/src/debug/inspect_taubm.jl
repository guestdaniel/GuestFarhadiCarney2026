using Helios

# Get WB params
Taumax = [0.0]
Taumin = [0.0]
ccall(
    (:Get_tauwb, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
    Cvoid,
    (
        Cint,
        Ptr{Cdouble},
        Cint,
        Cint,
        Ptr{Cdouble},
        Ptr{Cdouble},
    ),
    0,
    [1000.0],
    2,
    3,
    Taumax,
    Taumin,
)

# Get C1 params
bmTaumax = [0.0]
bmTaumin = [0.0]
ratiobm = [0.0]
ccall(
    (:Get_taubm, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
    Cvoid,
    (
        Cint,
        Ptr{Cdouble},
        Cint,
        Ptr{Cdouble},
        Ptr{Cdouble},
        Ptr{Cdouble},
        Ptr{Cdouble},
    ),
    0,
    [1000.0],
    2,
    Taumax,
    bmTaumax,
    bmTaumin,
    ratiobm,
)
bmTaubm = [1.0] .* (bmTaumax .- bmTaumin) .+ bmTaumin

# Get ???
TauWBMax = Taumin .+ 0.2 .* (Taumax .- Taumin)
TauWBMin = TauWBMax ./ Taumax .* Taumin
tauwb = TauWBMax .+ (bmTaubm .- bmTaumax) .* (TauWBMax .- TauWBMin) ./ (bmTaumax .- bmTaumin)

"Minimum τ    = $(Taumin[1] * 1000) ms"
"Maximum τ    = $(Taumax[1] * 1000) ms"
"Minimum C1 τ = $(bmTaumin[1] * 1000) ms"
"Maximum C1 τ = $(bmTaumax[1] * 1000) ms"
"C1 ratio     = $(ratiobm[1])"
"C1 τ = $(bmTaubm[1] * 1000) ms"
"Minumum WB τ = $(TauWBMin[1]*1000) ms"
"Maximum WB τ = $(TauWBMax[1]*1000) ms"
"WB τ = $(tauwb[1]*1000) ms"

