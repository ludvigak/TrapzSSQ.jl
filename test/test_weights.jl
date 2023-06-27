using Test
using TrapzSSQ
using FFTW

# Test that the weights computed using closed-form expressions
# correspond to the FFT of the exact integrals.
function weight_test(t0, m, N)
    _, k_vec = TrapzSSQ.fourier_coeffs(zeros(N))
    w_vec = TrapzSSQ.complex_pow_weights(t0, m, N)
    w_ref = fft(TrapzSSQ.complex_pow_integrals(t0, m, k_vec))/N
    @test isapprox(w_vec, w_ref, rtol=1e-10)
end

@testset "Complex power integrals" begin
    for N = 1000:1001 # odd and even
        @testset "N = $N" begin
            for t0 = [1-0.01im, 1+0.01im]
                @testset "t0 = $t0" begin
                    for m=1:4
                        @testset "m = $m" begin
                            weight_test(t0, m, N)
                        end
                    end
                end
            end
        end
    end
end
