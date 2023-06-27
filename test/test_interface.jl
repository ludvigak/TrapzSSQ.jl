include("src/geometry.jl")

function test_interface()
    N = 100
    h = 2*pi/N
    t_vec = collect(0:N-1)*h
    gamma(t) = cos(t) + 2im*sin(t) # ellipse
    gamma_vec = gamma.(t_vec)
    for t0 = [2.0 + 0.02im, 2.0 - 0.02im]
        z0 = gamma(t0)
        sigma = @. cos(t_vec)^2
        for m=[1, 2, "log"]
            @testset "t0 = $t0, m=$m" begin
                if m=="log"
                    integrand = @. sigma * log(abs(gamma(t_vec) - z0))
                    I_nossq = TrapzSSQ.eval_log_ssq(integrand, gamma_vec, t_vec, z0,
                                                    threshold=0)
                    w_nossq = TrapzSSQ.log_ssq_weights(gamma_vec, t_vec, z0,
                                                       threshold=0)
                    I_ssq = TrapzSSQ.eval_log_ssq(integrand, gamma_vec, t_vec, z0)
                    w_ssq = TrapzSSQ.log_ssq_weights(gamma_vec, t_vec, z0)
                else
                    integrand = @. sigma / (gamma(t_vec) - z0)^m                
                    I_nossq = TrapzSSQ.eval_pow_ssq(integrand, gamma_vec, t_vec, z0, m,
                                                    threshold=0)
                    w_nossq = TrapzSSQ.pow_ssq_weights(gamma_vec, t_vec, z0, m,
                                                       threshold=0)
                    I_ssq = TrapzSSQ.eval_pow_ssq(integrand, gamma_vec, t_vec, z0, m)
                    w_ssq = TrapzSSQ.pow_ssq_weights(gamma_vec, t_vec, z0, m)
                end
                T = h*sum(integrand)
                @testset "Zero threshold equals trapezoidal" begin
                    @test isapprox(I_nossq, T, rtol=1e-14)
                    @test isapprox(sum(integrand .* w_nossq), T, rtol=1e-14)
                end
                @testset "Interfaces match" begin
                    @test isapprox(sum(integrand .* w_ssq), I_ssq, rtol=1e-9)
                end
                @testset "Im(t0) large enough" begin
                    @test abs(I_ssq-T)/abs(I_ssq) > 1e-8
                end
                @testset "Integral is nonzero" begin
                    @test abs(I_ssq) > 0.1
                end
            end
        end
    end
end

@testset "Interface" begin
    test_interface()
end
;
