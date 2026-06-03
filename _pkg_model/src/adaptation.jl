export adapt_pla, adapt_pla_c, adapt_ea, adapt_ea_iir, adapt_ea_iir_parallel, adapt_ea_iir_parallel_sb, adapt_pla_clike, printparam

function adapt_pla_c(x, α, β; fs=100e3)
    y = zeros(size(x))
    I1 = [0.0]
    I2 = [0.0]
    for idx in eachindex(x)
        ccall(
            (:apply_powerlaw_adaptation, "C:\\Users\\dguest2\\cl_code\\Helios\\src\\model\\libgfc2023.so"), 
            Cvoid,
            (
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Cint,
                Cdouble,
                Cdouble,
                Cdouble,
                Cdouble,
                Cdouble,
                Ptr{Cdouble},
                Ptr{Cdouble},
            ),
            x,
            zeros(size(x)),
            I1,
            I2,
            idx-1,
            α,
            β,
            0.0,
            0.0,
            1/fs,
            y,
            zeros(size(x)),
        )
    end
    return y
end

function adapt_pla(x, α, β; fs=100e3)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * (1/fs) / ((n-j) * (1/fs) + β)
        end
    end
    return y, I
end

function adapt_pla_clike(x, α, β; fs=100e3)
    y = zeros(size(x))
    I = 0.0
    for n in eachindex(x)
        y[n] = max(0.0, x[n] - α * I)
        I = 0.0
        for j in 1:n
            j_c = j - 1
            n_c = n - 1
            I += y[j] * (1/fs) / ((n_c-j_c) * (1/fs) + β)
        end
    end
    return y
end

function adapt_ea_iir_parallel(x, α, τs; fs=100e3)
    y = zeros(length(x))
    I = zeros(length(x), length(τs))
    I_comb = zeros(length(x))
    ds = @. exp(-(1/fs)/τs)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I_comb[n-1])
        end

        # Compute I[n]
        for j in eachindex(ds)
            if n == 1
                I[n, j] = (1-ds[j]) * y[n]  # note: removed 1e-1 scaling factor?
            else
                I[n, j] = (1-ds[j]) * y[n] + ds[j] * I[n-1, j]
            end
        end
        I_comb[n] = sum(I[n, :])
    end
    return y, I_comb, I
end

function adapt_ea_iir_parallel_sb(x, α, τs; fs=100e3)
    y = zeros(length(x))
    I = zeros(length(τs))
    I_comb = 0.0
    ds = @. exp(-(1/fs)/τs)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I_comb)
        end

        # Compute I[n]
        for j in eachindex(ds)
            if n == 1
                I[j] = (1-ds[j]) * y[n]
            else
                I[j] = (1-ds[j]) * y[n] + ds[j] * I[j]
            end
        end
        I_comb = sum(I)
    end
    return y, I_comb, I
end


function adapt_ea(x, τₐ, τₑ; fs=100e3)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * 1/fs * exp( ((j-n) * (1/fs)) * (1/τₑ))
        end
    end
    return y, I
end

function adapt_ea_iir(x, τₐ, τ; fs=100e3)
    y = zeros(size(x))
    I = zeros(size(x))
    d = exp(-(1/fs)/τ)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        if n == 1
            I[n] = (1-d) * 1e-1 * y[n]
        else
            I[n] = (1-d) * 1e-1 * y[n] + d * I[n-1]
        end
    end
    return y, I
end

function printparam(θ)
    println("Best coefficient = $(round(θ[1]; digits=5))")
    println("Best short time constant = $(round(exp(θ[2]); digits=5)) ms")
    println("Best long time constant = $(round(exp(θ[3]); digits=5)) ms")
end

