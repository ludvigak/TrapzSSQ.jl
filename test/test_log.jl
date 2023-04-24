using Test
using TrapzSSQ
using QuadGK

# Numerical integration of the complex log is not straightforward,
# here we make do with some scalings and cherry-picked t0 to make
# it work.

function numeric_complex_log_integrals(t0, klist)
    N = length(klist);
    q = zeros(ComplexF64, N)
    if imag(t0) < 0
        s = -exp(1im*t0)
    else
        s = exp(1im*t0)
    end
    for i=1:N
        k = klist[i];
        f(t) = @. exp(1im*t)^k * (log(s*exp(1im*t)-s*exp(1im*t0)) - log(s))
        q[i], acc = quadgk(f, 0, 2*pi, atol=1e-10)
    end
    return q
end

@testset "Complex log integrals" begin
    for t0 = [1-0.1im, pi+0.1im]
        @testset "t0 = $t0" begin
            klist = -3:3
            q = TrapzSSQ.complex_log_integrals(t0, klist)
            qnum = numeric_complex_log_integrals(t0, klist)
            for (k, qk_test, qk_num) in zip(klist, q, qnum)
                @testset "k = $k" begin
                    if k != 0
                        @test isapprox(qk_test, qk_num, atol=1e-8)
                    else
                        @test isapprox(real(qk_test), real(qk_num), atol=1e-8)
                    end
                end
            end
        end
    end
end
;
