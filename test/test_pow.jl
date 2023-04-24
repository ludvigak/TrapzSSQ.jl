using Test
using TrapzSSQ
using QuadGK

function numeric_complex_pow_integrals(t0, m, klist)
    N = length(klist);
    p = zeros(ComplexF64, N)       
    for i=1:N
        k = klist[i];
        f(t) = @. exp(1im*k*t) / (exp(1im*t)-exp(1im*t0))^m
        p[i], _ = quadgk(f, 0, 2*pi)
    end
    return p
end

@testset "Complex power integrals" begin
    for t0 = [1-0.1im, 1+0.1im]
        @testset "t0 = $t0" begin
            for m=1:2
                @testset "m = $m" begin
                    kmax = m
                    klist = -kmax:kmax       
                    p = TrapzSSQ.complex_pow_integrals(t0, m, klist)
                    pnum = numeric_complex_pow_integrals(t0, m, klist)
                    for (k, ptest, pref) in zip(klist, p, pnum)
                        @testset "k = $k" begin
                            @test isapprox(ptest, pref, atol=1e-10)
                        end
                    end
                end
            end
        end
    end
end        
;
