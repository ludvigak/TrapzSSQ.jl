struct CurveDefinition
    gamma   :: Function
    gammap  :: Function
    gammapp :: Function
end

struct CurveDiscretization
    N       :: Integer
    h       :: Real
    t       :: Vector
    gamma   :: Vector
    gammap  :: Vector
    gammapp :: Vector
end

function starfish(amplitude, n_arms, radius=1.0)
    gamma(t) = radius*(1 + amplitude*cos(n_arms*t))*exp(1im*t)
    gammap(t) = radius*((1 + amplitude*cos(n_arms*t))*exp(1im*t)*1im +
         (-n_arms*amplitude*sin(n_arms*t)).*exp(1im*t))
    gammapp(t) = -radius*exp(1im*t)*(amplitude*cos(n_arms*t) + amplitude*n_arms^2*cos(n_arms*t) + 1 + amplitude*n_arms*sin(n_arms*t)*2im)
    return CurveDefinition(gamma, gammap, gammapp)
end

function discretize(curve::CurveDefinition, N::Integer)
    h = 2*pi/N
    t_vec = collect(0:N-1)*h
    return CurveDiscretization(
        N,
        h,
        t_vec,
        curve.gamma.(t_vec),
        curve.gammap.(t_vec),
        curve.gammapp.(t_vec),
    )
end

function reference_problem(curve :: CurveDefinition, t0, z0, m)
    @assert abs((z0-curve.gamma(t0))) < 2*eps() # Sanity check
    # Define kernel and integrand
    if m=="log"
        integrand = (sigma, gamma, gammap, z0) -> sigma * abs(gammap) * log(abs(gamma-z0))
    else
        integrand = (sigma, gamma, gammap, z0) -> sigma * gammap / (gamma-z0)^m
    end
    # Define density and reference solution
    I_ana = nothing # Default is we don't know analytical solution
    if m=="log"
        sigma = t -> real(curve.gamma(t)) * imag(curve.gamma(t))
    else
        if imag(t0) < 0
            # z0 outside curve
            field  = z -> 1/z
            I_ana = 1/(-z0)^m
        else
            # z0 inside curve
            field   = z -> z^3 + z
            fieldp  = z -> 3*z^2 + 1
            fieldpp = z -> 6*z
            if m==1
                I_ana = field(z0)
            elseif m==2
                I_ana = fieldp(z0)
            elseif m==3
                I_ana = fieldpp(z0)/2
            elseif m==4
                I_ana = 1
            else
                I_ana = 0
            end
        end
        sigma = t -> field(curve.gamma(t)) / (2im*pi)
    end
    integrand, sigma, I_ana
end