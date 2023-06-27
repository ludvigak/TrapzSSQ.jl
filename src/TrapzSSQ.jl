module TrapzSSQ

using Printf

include("fourier.jl")
include("weights.jl")

function find_root(z0, gamma_vec, t_vec, expansion::FourierSeries)
    minidx = argmin(@. abs(gamma_vec-z0))
    tinit = t_vec[minidx]
    t0 = complex_fourier_newton(expansion, z0, tinit)
    return t0
end

function eval_pow_ssq(integrand_vec::Vector, gamma_vec::Vector, t_vec::Vector, z0::Complex, m::Integer; gamma_hat=FourierSeries(gamma_vec), threshold=Inf)
    t0 = find_root(z0, gamma_vec, t_vec, gamma_hat)
    if abs(imag(t0)) > threshold
        # Return trapezoidal sum
        return sum(integrand_vec)*2*pi/length(integrand_vec)
    end
    f_vec = @. integrand_vec * (exp(1im*t_vec) - exp(1im*t0))^m
    f_hat = FourierSeries(f_vec)
    p_vec = TrapzSSQ.complex_pow_integrals(t0, m, f_hat.k)
    I_ssq = sum( p_vec .* f_hat.coeffs )
    return I_ssq
end

function eval_log_ssq(integrand_vec::Vector{<:Real}, gamma_vec::Vector, t_vec::Vector, z0::Complex; gamma_hat=FourierSeries(gamma_vec), threshold=Inf)
    # Note that integrand must be real
    N = length(gamma_vec)
    h = 2*pi/N
    t0 = find_root(z0, gamma_vec, t_vec, gamma_hat)
    if abs(imag(t0)) > threshold
        # Return trapezoidal sum
        return h * sum(integrand_vec)
    end
    log_vec = @. log(abs(gamma_vec-z0))
    mapped_log_vec = @. log(abs(exp(1im*t_vec) - exp(1im*t0)))
    f_vec = @. integrand_vec / log_vec
    rem_vec = @. f_vec * (log_vec - mapped_log_vec)
    f_hat = FourierSeries(f_vec)
    q_vec = TrapzSSQ.complex_log_integrals(t0, f_hat.k)
    I_ssq = h * sum(rem_vec) + real(sum( q_vec .* f_hat.coeffs ))
    return I_ssq
end

function pow_ssq_weights(gamma_vec::Vector, t_vec::Vector, z0::Complex, m::Integer;
                         gamma_hat=FourierSeries(gamma_vec), threshold=Inf)
    t0 = find_root(z0, gamma_vec, t_vec, gamma_hat)
    N = length(gamma_vec)
    if abs(imag(t0)) > threshold
        # Return trapezoidal sum
        return 2*pi/N * ones(N)
    end
    w_vec = complex_pow_weights(t0, m, N)
    ssq_weights = @. w_vec * (exp(1im*t_vec) - exp(1im*t0))^m
    return ssq_weights
end

function log_ssq_weights(gamma_vec::Vector, t_vec::Vector, z0::Complex;
                         gamma_hat=FourierSeries(gamma_vec), threshold=Inf)
    N = length(gamma_vec)
    h = 2*pi/N
    t0 = find_root(z0, gamma_vec, t_vec, gamma_hat)
    if abs(imag(t0)) > threshold
        # Return trapezoidal rule weights
        return h * ones(N)
    end
    inv_log_vec = @. 1 / log(abs(gamma_vec-z0))
    mapped_log_vec = @. log(abs(exp(1im*t_vec) - exp(1im*t0)))
    w_vec = real(fft( TrapzSSQ.complex_log_integrals(t0, gamma_hat.k) ))/N
    ssq_weights = @. h*(1 - mapped_log_vec*inv_log_vec) + w_vec*inv_log_vec
    return ssq_weights
end

function complex_pow_integrals(t0, m, klist)
    N = length(klist)
    p = zeros(ComplexF64, N)
    mfact = factorial(m-1)
    for i=1:N
        k = klist[i];
        if imag(t0)>0 && k >= m
            sgn = 1;
        elseif imag(t0)<0 && k<=0
            sgn = -1;
        else
            sgn = 0;
        end
        c = prod(k .- (1:m-1));   
        val = 2*pi*c/mfact*exp(1im*(k-m)*t0)
        p[i] = sgn*val
    end
    return p
end

function complex_log_integrals(t0, klist)
    N = length(klist)    
    q = zeros(ComplexF64, N);
    # Non-vectorized
    for (i, k) in enumerate(klist)
        qk = 0
        if imag(t0)<0
            if k<0
                qk = 2*pi/k*exp(1im*k*t0)
            elseif k==0
                qk = -2*pi*imag(t0)
            else # k>0
                qk = 0
            end        
        elseif imag(t0)>0
            if k<0
                qk = 2*pi/k
            elseif k==0
                qk = 0
            else # k>0
                qk = 2*pi/k*(1-exp(1im*k*t0))
            end
        end
        q[i] = qk
    end
    return q
end

end
