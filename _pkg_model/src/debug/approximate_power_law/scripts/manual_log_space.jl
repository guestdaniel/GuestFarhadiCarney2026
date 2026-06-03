tau_short = 5e-3
tau_long = 1e4
n = 125

delta = (log(tau_long) - log(tau_short))/(n-1)

j = collect(LogRange(tau_short, tau_long, n))
c = zeros(n)
for ii in eachindex(c)
    ii_c = ii - 1
    c[ii] = exp(log(tau_short) + delta * ii_c)
end