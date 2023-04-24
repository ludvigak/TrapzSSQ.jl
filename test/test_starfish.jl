using Test
using TrapzSSQ
using QuadGK

include("src/geometry.jl")

function starfish_error(t0ref, m, N)
    # Define geometry
    curve = starfish(0.3, 5)
    z0 = curve.gamma(t0ref)
    # Define kernel and integrand
    integrand, sigma, I_ana = reference_problem(curve, t0ref, z0, m)
    ###### Compute
    # Discretize
    discretization = discretize(curve, N)
    sigma_vec = sigma.(discretization.t)
    # Compute trapezoidal quadrature
    int_vec = integrand.(sigma_vec, discretization.gamma, discretization.gammap, z0)
    I_trapz = discretization.h * sum(int_vec)
    ### Begin Singularity Swap Quadrature (SSQ) 
    if m=="log"
        I_ssq = TrapzSSQ.eval_log_ssq(int_vec, discretization.gamma, discretization.t, z0)
    else
        I_ssq = TrapzSSQ.eval_pow_ssq(int_vec, discretization.gamma, discretization.t, z0, m)
    end
    ### End SSQ
    if !isnothing(I_ana)
        Erel_ssq = abs(I_ssq - I_ana) / abs(I_ana)
    else
        # Compute reference numerically
        analytic_integrand = t -> integrand(sigma(t), curve.gamma(t), curve.gammap(t), z0)
        I_ref, err_ref = quadgk(analytic_integrand, 0, 2*pi, rtol=1e-13)
        Erel_ssq = abs(I_ssq - I_ref) / abs(I_ref)
    end
    return Erel_ssq
end

@testset "Starfish reference" begin
    N = 401 # odd
    m_list =   ["log", 1,     2,     3    , 4,  ]
    tol_list = [1e-13, 1e-14, 1e-12, 1e-10, 1e-9]
    for (m, tol) = zip(m_list, tol_list)
        @testset "m=$m" begin
            for t0ref = [1.0 + 0.01im, 1.0 - 0.01im]
                @testset "t0=$t0ref" begin
                    @test starfish_error(t0ref, m, N) < tol
                end
            end
        end
    end
end