# Solving 1D diffusion equation to steady state with central finite differences
using CairoMakie
# Define Function
@views function steadyDiffusion_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    ρ       = (lx/(dc*2π))^2        # Numerical density
    # Numerics
    ncx     = 200                   # Number of cells in x
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                  # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)      # Convergence check frequency
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx/sqrt(1/ρ) #dx^2/dc/2.1           # Time step size
    nt      = 5ncx #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    fig1    = Figure()                 # Plotting
    ax1     = Axis(fig1[1, 1])
    ax2     = Axis(fig1[2, 1], yscale = log10)
    lines!(ax1, xc, C, color = :blue)
    lines!(ax1, xc, C_ini, color = :orange)
    lines!(ax2, iter_evo, err_evo, color = :blue)
    display(fig1) 
    # Iteration loop
    while err >= ϵtol && iter <= maxiter
        # Computation
        qx         .-=   dt./(ρ.*dc + dt).*(qx + dc.*diff(C)./dx)                  # Ludo's formulation 
        C[2:end-1] .-=   dt.* diff(qx)./dx
        # Visualisation
        if iter % ncheck == 0
            err = maximum(abs.(diff(dc.*diff(C)./dx)./dx))
            push!(iter_evo, iter/ncx); push!(err_evo, err)
            sleep(0.1)
            empty!(ax1)
            empty!(ax2)
            lines!(ax1, xc, C, color = :blue)
            lines!(ax1, xc, C_ini, color = :orange)
            lines!(ax2, iter_evo, err_evo, color = :blue)
            display(fig1)
        end
        iter += 1
    end        
end

# Call function
steadyDiffusion_implicit_1D()