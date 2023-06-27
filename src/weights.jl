#
# Routines for computing modified quadrature weights directly
#

using Symbolics

# Dictionary for keeping generated functions
Sm_func_table = Dict{Integer, Tuple{Function, Function}}()

function generate_Sm_func(m)
    # Symbolic generation of Sm, following paper
    @variables r, Kmin, Kmax
    Dr = Differential(r)
    sum_impos = (r^Kmax - r^(m-1))/(r-1)
    sum_imneg = (r^(Kmin-1)-1)/(r-1)
    if m > 1
        Sm_impos = r^m * expand_derivatives((Dr^(m-1))(sum_impos))
        Sm_imneg = r^m * expand_derivatives((Dr^(m-1))(sum_imneg))
    else
        Sm_impos = r * sum_impos
        Sm_imneg = r * sum_imneg
    end
    Sm_impos_expr = build_function(Sm_impos, r, Kmin, Kmax, expression=Val{false})
    Sm_imneg_expr = build_function(Sm_imneg, r, Kmin, Kmax, expression=Val{false})
    # Why these?
    Base.remove_linenums!(Sm_impos_expr)    
    Base.remove_linenums!(Sm_imneg_expr)
    # Return both functions
    return (Sm_impos_expr, Sm_imneg_expr)
end

function get_Sm_func(t0, m::Integer)
    # Table lookup with lazy generation
    if !haskey(Sm_func_table, m)
        @info "Generating exact weights for m=$m"
        Sm_func_table[m] = generate_Sm_func(m)
    end
    (Sm_func_impos, Sm_func_imneg) = Sm_func_table[m]
    if imag(t0)>0
        return Sm_func_impos
    else
        return Sm_func_imneg
    end
end

function complex_pow_weights(t0, m, N)
    w_vec = zeros(ComplexF64, N)
    Kmin, Kmax = k_range(N)
    # Vectorized implementation
    i=0:N-1
    ti = 2*pi*i/N
    x = @. t0-ti
    r = @. exp(1im*x)
    # Hardcoded expressions for m=1 and m=2, use generated for higher m
    if m==1
        if imag(t0)>0
            Sm = @. r*(r^Kmax - 1) / (r - 1)
        else
            Sm = @. r*(r^(Kmin-1)-1) / (r-1)
        end
    elseif m==2
        if imag(t0)>0
            Sm = @. r^2*((Kmax*r^(Kmax - 1) - 1)/(r - 1) + (r - r^Kmax)/(r - 1)^2)
        else
            Sm = @. (r^Kmin - Kmin*r^Kmin + r^2 - 2*r*r^Kmin + Kmin*r*r^Kmin)/(r - 1)^2
        end
    else
        Sm_func = get_Sm_func(t0, m)
        Sm = Sm_func.(r, Kmin, Kmax)    
    end
    w_vec = 2*pi/N*exp(-1im*m*t0)/factorial(m-1)*Sm        
    return w_vec
end


