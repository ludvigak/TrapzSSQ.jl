# TrapzSSQ: Singularity swap quadrature for the trapezoidal rule

This repository accompanies the paper "Singularity swap quadrature for nearly singular line integrals on closed curves in two dimensions" by L. af Klinteberg (arXiv).

The code is a Julia implementation of the quadrature method introduced in the paper. All plots reported in the paper can be generated using the scripts in `scripts/`.

### Example usage

Start Julia from the repo root with `julia --project=.`

Install alla dependencies using
```julia
using Pkg; Pkg.instantiate() # First time setup
Pkg.test("TrapzSSQ")         # Test that it works
```

Run the demo scripts:
```julia
include("scripts/demo_laplace.jl")
include("scripts/show_decay.jl")
include("scripts/show_convergence.jl")
```
This will output the plot files in the current directory as png and pdf files.
