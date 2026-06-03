using CairoMakie
using AuditorySignalUtils

function fc_to_d(fc)
    exp(-2π * fc/100e3)
end

function fc_to_d_old(fc)
    exp(-2π * fc/50e3)
end

function d_to_fc(d)
    -log(d)/(2π)*100e3
end

function d_to_fc_old(d)
    -log(d)/(2π)*50e3
end

# OKAY
fig = Figure()
d_afagh = 1.0 - 3.9998e-5
fc = LogRange(1e-2, 2e0, 1000)
ax = Axis(fig[1, 1], xscale=log10)
scatter!(ax, [0.2], [fc_to_d(0.2)]; color=:black)
scatter!(ax, [0.2], [fc_to_d_old(0.2)]; color=:black, marker=:rect)
scatter!(ax, [d_to_fc(d_afagh)], [d_afagh]; color=:red)
scatter!(ax, [d_to_fc_old(d_afagh)], [d_afagh]; color=:red, marker=:rect)
lines!(ax, fc, fc_to_d.(fc))
ax.xticks = vcat((1:10) .* (10.0^-1), (1:10) .* (10.0^0))
fig