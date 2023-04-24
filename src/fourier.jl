using FFTW

struct FourierSeries
    coeffs :: Vector{ComplexF64}
    k      :: Vector{Float64}
    FourierSeries(values) = new(fourier_coeffs(values)...)
end

function fourier_coeffs(fj)
    fj = fj
    N = length(fj)
    k = ifftshift(get_k_vec(N, 2*pi))
    c = fft(fj)/N
    return c, k
end

function fourier_magnitude(fj)
    cr, k = fourier_coeffs(real(fj))
    ci, k = fourier_coeffs(imag(fj))
    m = @. sqrt(abs(cr)^2 + abs(ci)^2)
    return fftshift(m), fftshift(k)
end
    
function get_k_vec(M,L)
    if mod(M,2)==0
        MM = M/2
        k = (2*pi/L)*collect(-MM:(MM-1));
    elseif mod(M-1,2)==0
        MM = (M-1)/2
        k = (2*pi/L)*collect(-MM:MM)
    else 
        error("k-vectors not computed");
    end
    return k
end

function eval_expansion(expansion::FourierSeries, t)
    return eval_expansion(expansion.coeffs, expansion.k, t)
end

function eval_expansion(c::Vector{ComplexF64}, k::Vector{Float64}, t)
    e = exp.(1im*k*t)
    f = sum(c .* e)
    fp = sum(1im * k .* c .* e)
    return f, fp
end

function complex_fourier_newton(expansion::FourierSeries, z::ComplexF64, tinit; tol=1e-13, maxiter=50)
    t = ComplexF64(tinit)
    dt = ComplexF64(0)
    for iter=1:maxiter
        gamma, gammap = eval_expansion(expansion, t)
        F = gamma - z
        Fp = gammap
        dt = -F/Fp
        t = t + dt
        if abs(dt) < tol
            return t
        end
    end
    # Did not converge
    @warn(@sprintf("Newton failed, |dt|=%.3e, t=%g+%gi", abs(dt), real(t), imag(t)))
    return t
end
