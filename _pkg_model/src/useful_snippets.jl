# Useful code snippets for diagnosing things in the model

# 1) Plot Greenwood function (CF to mm in cat);
# From line ~450 in `model.c`
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3]); 
lines!(ax, cf, 11.9 .* log10.(0.80 .+ cf ./ 456.0)); 
ax.xlabel = "CF (Hz)"; ax.ylabel = "Distance (mm)";
fig


# 2) Plot CF versus control-path CF
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3], yscale=log10, yticks=[0.1e3, 1e3, 10e3])
bmplace = 11.9 .* log10.(0.80 .+ cf ./ 456.0)
cf_shifted = 456.0 .* (10 .^ ((bmplace .+ 1.2) ./ 11.9) .- 0.80)
lines!(ax, cf, cf_shifted)
ax.xlabel = "CF (Hz)"; ax.ylabel = "CF (control path, Hz)";
ylims!(ax, 0.1e3, 40e3)
xlims!(ax, 0.1e3, 40e3)
fig

# 3) Plot gain versus CF
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3])
gain = @. 52.0/2.0 * (tanh(2.2*log10(cf/0.6e3)+0.15) + 1.0)
lines!(ax, cf, gain; color=:black)
lines!(ax, cf, max.(gain, 15.0); color=:black, linestyle=:dash)
xlims!(ax, 0.1e3, 40e3)
ax.xlabel = "CF (Hz)"; ax.ylabel = "Gain (dB)";
fig

# 4) Plot Q10 versus CF
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3])
Q10 = @. 10.0 ^ (0.4708*log10(cf/1e3) + 0.4664)
lines!(ax, cf, Q10; color=:black)
xlims!(ax, 0.1e3, 40e3)
ylims!(ax, 0.0, 20.0)
ax.xlabel = "CF (Hz)"; ax.ylabel = "Q10";
fig

# 5) Plot BW versus CF
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3])
Q10 = @. 10.0 ^ (0.4708*log10(cf/1e3) + 0.4664)
lines!(ax, cf, cf ./ Q10; color=:black)
xlims!(ax, 0.1e3, 40e3)
ax.xlabel = "CF (Hz)"; ax.ylabel = "Bandwidth (Hz)";
fig

# 6) Plot taumax vs taumin
cf = LogRange(0.1e3, 40e3, 1000)
fig = Figure(); 
ax = Axis(fig[1, 1]; xscale=log10, xticks=[0.1e3, 1e3, 10e3], yscale=log10, yticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
Q10 = @. 10.0 ^ (0.4708*log10(cf/1e3) + 0.4664)
bw = cf ./ Q10
gain = @. 52.0/2.0 * (tanh(2.2*log10(cf/0.6e3)+0.15) + 1.0)
gain = max.(gain, 15.0)
τ_max = @. 2.0 / (2π * bw)
τ_min = @. τ_max * 10.0 ^ (-gain/(20.0*3.0))
lines!(ax, cf, τ_min .* 1e3; color=:black)
lines!(ax, cf, τ_max .* 1e3; color=:black, linestyle=:dash)
xlims!(ax, 0.1e3, 40e3)
ax.xlabel = "CF (Hz)"; ax.ylabel = "τ WB (ms)"
vlines!(ax, [600.0]; color=:gray, linestyle=:dash)
fig
