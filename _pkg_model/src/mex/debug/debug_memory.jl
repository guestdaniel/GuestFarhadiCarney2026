using Helios
using CairoMakie

function meminfo_julia()
    Base.gc_live_bytes()/2^20
end

mems = map(1:25) do rep
    sim_gfc2023(zeros(50000), [1000.0])
    meminfo_julia()
end

lines(mems)
