using TrapzSSQ
using Plots
using FFTW
using Printf
using QuadGK
using LinearAlgebra
using LaTeXStrings

include("../test/src/geometry.jl")

###### Setup

function show_decay(t0ref)
    # Define geometry
    curve = starfish(0.3, 5, 1)
    @show t0ref
    z0 = curve.gamma(t0ref)
    @show z0

    interior = (imag(t0ref) > 0)
    @show interior

    # Define discretization size
    N = 401 # odd

    # Define kernel and integrand
    m = 1
    integrand, sigma, I_ana = reference_problem(curve, t0ref, z0, m)
    ###### Compute

    # Discretize
    discretization = discretize(curve, N)
    sigma_vec = sigma.(discretization.t)

    # Compute quadrature
    int_vec = integrand.(sigma_vec, discretization.gamma, discretization.gammap, z0)
    I_trapz = discretization.h * sum(int_vec)
    I_ssq = TrapzSSQ.eval_pow_ssq(int_vec, discretization.gamma,discretization.t, z0, m)

    ###### Report
    @show m
    @show I_ssq
    if !isnothing(I_ana)
        @show abs(I_trapz - I_ana)
        @show abs(I_ssq - I_ana)
    end

    # Compute some of the SSQ quantities so that we can show them
    # Compute Fourier expansion of discretization
    gamma_hat = TrapzSSQ.FourierSeries(discretization.gamma)
    t0 = TrapzSSQ.find_root(z0, discretization.gamma, discretization.t, gamma_hat)
    @show t0
    f_vec = @. int_vec * (exp(1im*discretization.t) - exp(1im*t0))^m
    f_mag, k_plot = TrapzSSQ.fourier_magnitude(f_vec)
    int_mag, k_plot = TrapzSSQ.fourier_magnitude(int_vec)
    plot(xlabel="wavenumber, "*L"k", ylabel="Fourier coefficient magnitude",
        yaxis=:log10, reuse=false,
        xticks=-200:100:200,
        yticks=[1e-15, 1e-10, 1e-5, 1e0])
    plot!(k_plot, f_mag, label=L"|F[f](k)|")
    plot!(k_plot, int_mag, label=L"|F[\sigma (\tau-z)^{-1}](k)|")
    rho = exp(abs(imag(t0)))
    plot!(k_plot, 0.25*rho.^(-abs.(k_plot)), style=:dash, 
        label=L"O(\rho^{-|k|})")
end

gr()
plot_font = "Computer Modern"
default(fontfamily=plot_font, 
    legendfontsize=12,
    labelfontsize=12,
    linewidth=2, gridwidth=2, framestyle=:box, label=nothing, grid=true
    )
default(legend=:topright)
show_decay(1.0 + 0.05im)
savefig("decay1.pdf")
println("Saved decay1.pdf")
default(legend=:right)
show_decay(1.0 + 0.005im)
savefig("decay2.pdf")
println("Saved decay2.pdf")
