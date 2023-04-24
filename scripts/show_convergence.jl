using TrapzSSQ
using LinearAlgebra
using Plots
using LaTeXStrings
using QuadGK

include("../test/src/geometry.jl")

function run_convergence(N_vec, t0_vec, m; runlabel="")
    curve = starfish(0.3, 5, 1)
    max_errors_tpz = zeros(length(N_vec))
    max_errors_ssq = zeros(length(N_vec))
    I_ref = zeros(ComplexF64, length(t0_vec))
    meandist = 0
    for (i, t0) in enumerate(t0_vec)
        z0 = curve.gamma(t0)
        znb = curve.gamma(real(t0))
        meandist += abs(z0-znb) / length(t0_vec)
        integrand, sigma, I_ana = reference_problem(curve, t0, z0, m)
        if isnothing(I_ana)
            if i==1
                println(" * Using numerical reference")
            end 
            # Compute reference
            analytic_integrand = t -> integrand(
                sigma(t), curve.gamma(t), curve.gammap(t), z0)
            I_ref[i], err_ref = quadgk(analytic_integrand, 0, 2*pi, rtol=1e-15)
        else
            I_ref[i] = I_ana
        end
    end
    @show meandist
    I_norm = maximum(abs.(I_ref))
    for (N_idx, N) in enumerate(N_vec)
        discretization = discretize(curve, N)   
        gamma_hat = TrapzSSQ.FourierSeries(discretization.gamma)
        
        errors_tpz = zeros(length(t0_vec))
        errors_ssq = zeros(length(t0_vec))
        for (t0_idx, t0) in enumerate(t0_vec)
            z0 = curve.gamma(t0)
            integrand, sigma, I_ana = reference_problem(curve, t0, z0, m)
            sigma_vec = sigma.(discretization.t)
            int_vec = integrand.(sigma_vec, discretization.gamma, discretization.gammap, z0)
            I_tpz = discretization.h * sum(int_vec)
            if m=="log"
                I_ssq = TrapzSSQ.eval_log_ssq(
                    int_vec, discretization.gamma,discretization.t, z0, gamma_hat=gamma_hat)
            else
                I_ssq = TrapzSSQ.eval_pow_ssq(
                    int_vec, discretization.gamma,discretization.t, z0, m, gamma_hat=gamma_hat)
            end
            errors_tpz[t0_idx] = abs(I_tpz - I_ref[t0_idx])
            errors_ssq[t0_idx] = abs(I_ssq - I_ref[t0_idx])
        end
        max_errors_tpz[N_idx] = maximum(errors_tpz) / I_norm + eps()
        max_errors_ssq[N_idx] = maximum(errors_ssq) / I_norm + eps()
    end

    Emax = max(1, maximum(max_errors_tpz), maximum(max_errors_ssq))
    plot!(yaxis=:log10, ylim=(1e-16, Emax), yticks=[1e-15, 1e-10, 1e-5, 1],
        ylabel="maximum relative error", xlabel="N", title="Kernel: " * (m=="log" ? "log" : "m=$m"))
    plot!(N_vec, max_errors_tpz, label=runlabel*"Trapezoidal",
        markershape=:square, markersize=3, markerstrokewidth=0)
    plot!(N_vec, max_errors_ssq, label=runlabel*"SSQ", 
        markershape=:circle, markersize=3, markerstrokewidth=0)
end

function main()
    for sign=[-1, 1]
        for m=[1, 2, 3, "log"]
            println("================== $m ================")
            N_vec = 25:25:1000
            t0list = LinRange(0, 2*pi, 100)

            gr()
            plot_font = "Computer Modern"
            default(fontfamily=plot_font, 
                legendfontsize=12,
                labelfontsize=12,
                linewidth=2, gridwidth=2, framestyle=:box, label=nothing, grid=true,
                legend=:topright
                )
            plot()

            # interior if d > 0
            d = sign*0.01
            run_convergence(N_vec, t0list .+ d*1im, m, runlabel="d=$d, ")

            d = sign*0.02
            run_convergence(N_vec, t0list .+ d*1im, m, runlabel="d=$d, ")

            d = sign*0.04
            run_convergence(N_vec, t0list .+ d*1im, m, runlabel="d=$d, ")
            if true # Save figues?
                filename = "convergence_$m" *
                    (sign > 0 ? "_interior" : "_exterior") * ".pdf"
                savefig(filename)
                println("Saved $filename")
            else
                gui()
            end
        end
    end
end

main()