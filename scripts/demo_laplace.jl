using TrapzSSQ
using Plots
using LinearAlgebra
using LaTeXStrings

include("../test/src/geometry.jl")

@inline
function dlp_integrand(grid::CurveDiscretization, sigma::Vector{Float64}, target::ComplexF64) :: Vector{ComplexF64}
    return @. sigma * grid.gammap/(grid.gamma-target) / (2*pi)
end

function dlp_eval(grid::CurveDiscretization, sigma::Vector{Float64}, target::ComplexF64) :: Float64
    return grid.h * sum(imag(dlp_integrand(grid, sigma, target)))
end

function dlp_matrix(grid)
    A = zeros(Float64, grid.N, grid.N)
    for i=1:grid.N
        for j=1:grid.N
            A[i,j] = grid.h*imag( grid.gammap[j]/(grid.gamma[j]-grid.gamma[i]) )/(2*pi);
        end
        A[i,i] = 1/2 + grid.h*imag(grid.gammapp[i]/(2*grid.gammap[i]))/(2*pi);
    end    
    return A
end

function main()
    curve = starfish(0.3, 5, 1)
    N = 400
    threshold=0.1
    discretization = discretize(curve, N)
    A = dlp_matrix(discretization)
    zsing = 3+3im
    fexact(z) = log(abs(z-zsing))
    rhs = fexact.(discretization.gamma)
    sigma_vec = A\rhs
    # Evaluation grid
    Mx = 400
    My = 400
    xgrid = collect(LinRange(-1.299, 1.299, Mx))
    ygrid = collect(LinRange(-1.299, 1.299, My))
    DLP = zeros(Mx, My)
    REF = zeros(Mx, My)
    SSQ = zeros(Mx, My)
    L = discretization.h * sum(abs.(discretization.gammap))
    ds = L/N
    one_vec = ones(N)
    gamma_hat = TrapzSSQ.FourierSeries(discretization.gamma)
    for i=1:Mx
        for j=1:My
            z = xgrid[i] + ygrid[j]*1im
            REF[j,i] = fexact(z)
            # Figure out if interior point
            dist = minimum(abs.(discretization.gamma .- z))
            if dist > 2*ds
                if dlp_eval(discretization, one_vec, z) < 0.5
                    DLP[j,i] = NaN
                    SSQ[j,i] = NaN
                    continue
                end
            else
                marker_int = dlp_integrand(discretization, one_vec, z)
                marker = imag(TrapzSSQ.eval_pow_ssq(marker_int, discretization.gamma, discretization.t, z, 1, gamma_hat=gamma_hat, threshold=threshold))
                if marker < 0.5
                    DLP[j,i] = NaN
                    SSQ[j,i] = NaN
                    continue
                end
            end
            # Compute
            DLP[j,i] = dlp_eval(discretization, sigma_vec, z)
            if dist > 8*ds
                SSQ[j,i] = DLP[j,i]
            else
                int_vec = dlp_integrand(discretization, sigma_vec, z)
                SSQ[j,i] = imag(
                    TrapzSSQ.eval_pow_ssq(
                        int_vec, discretization.gamma, 
                        discretization.t, z, 1, 
                        gamma_hat=gamma_hat,
                        threshold=threshold
                        )
                    )
            end
        end
    end
    ERR_TPZ = @. abs(DLP-REF) + eps(REF)
    ERR_SSQ = @. abs(SSQ-REF) + eps(REF)

    function plot_and_save_error(E, filename)
        heatmap(xgrid, ygrid, log10.(E), aspect_ratio=1,
            clim=(-16,-0),
            xticks=[-1.3, 0, 1.3], yticks=[-1.3, 0, 1.3], 
            c=:diverging_bwr_20_95_c54_n256,
            colorbar_title=L"\log_{10} E",
            colorbar_ticks=[-15, -10, -5, 0],
            size=(450,400),
            dpi=300,
            interpolation=false
            )
        tplot = LinRange(0, 2*pi, 1000)
        cplot = curve.gamma.(tplot)
        plot!(real(cplot), imag(cplot), 
            width=2, color=:black)
        savefig(filename)
        println("Write $filename")
    end
    plot_and_save_error(ERR_TPZ, "laplace_tpz.png")
    plot_and_save_error(ERR_SSQ, "laplace_ssq.png")

    mask = @. !isnan(SSQ)
    maxrelerr_ssq = norm(ERR_SSQ[mask], Inf)  / norm(REF[mask], Inf)
    maxrelerr_tpz = norm(ERR_TPZ[mask], Inf)  / norm(REF[mask], Inf)

    @show maxrelerr_ssq
    @show maxrelerr_tpz
end

# pyplot() # Need PyPlot to get the colorbar labels right, otherwise GR is fine
gr()
#plot_font = "Computer Modern"
default(
    #fontfamily=plot_font, 
    legendfontsize=12,
    labelfontsize=12,
    linewidth=2, gridwidth=2, framestyle=:box, label=nothing, grid=false
    )

@time main()
#gui()