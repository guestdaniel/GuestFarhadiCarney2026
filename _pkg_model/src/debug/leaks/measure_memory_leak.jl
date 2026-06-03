## measure_memory_leak.m
# This script is intended to repeatedly run the efferent model, either by
# single for loop or by parfor, and to look for evidence of a memory leak.
# Current memory usage will be recorded and logged during runs so that 
# memory usage can be reconstructed afterward even if a crash occurs.

using Dates
using Printf
using Helios
using CairoMakie

function meminfo_julia()
  # @printf "GC total:  %9.3f MiB\n" Base.gc_total_bytes(Base.gc_num())/2^20
  # Total bytes (above) usually underreports, thus I suggest using live bytes (below)
  #Base.gc_live_bytes()/2^20
  #Base.jit_total_bytes()/2^20
  Sys.maxrss()/2^20
end

# Loop over series of empty efferent model runs
# Plot
fig = Figure()
ax = Axis(fig[1, 1])

# Run simulations
data = map(1:500) do _
	sim_gfc2023(zeros(50000), 1000.0)
    meminfo_julia()
end
lines!(ax, 1:length(data), data)
ax.xlabel = "Repeat #"
ax.ylabel = "Memory usage (MB)"
fig
