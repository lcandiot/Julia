# Solving 1D diffusion reaction equation to steady state with central finite differences
using Pkg, CairoMakie
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
# Define Function
@views function diffusionReaction_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    ρ       = (lx/(dc*2π))^2        # Numerical density
    Ceq     = 0.1                   # Equilibrium concentration
    da      = 10.0                  # Damkoehler no.
    ξ       = lx^2/dc/da
    # Numerics
    ncx     = 200                   # Number of cells in x
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                  # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)      # Convergence check frequency
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dτ      = dx/sqrt(1/ρ) #dx^2/dc/2.1           # Time step size
    nt      = 5ncx #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    fig1    = Figure()                 # Plotting
    ax1     = Axis(fig1[1, 1], xlabel=L"\textit{x}[m]", ylabel=L"\textit{C}[mol]", limits=(0, lx, 0, 2), title="1D Steady-state Diffusion (impl.)")
    ax2     = Axis(fig1[2, 1], yscale = log10, xlabel=L"\textit{iter/ncx}", ylabel=L"||\textit{r}||_{∞}")
    lines!(ax1, xc, C, color = :blue)
    lines!(ax1, xc, C_ini, color = :orange)
    lines!(ax2, iter_evo, err_evo, color = :blue)
    display(fig1) 
    # Iteration loop
    while err >= ϵtol && iter <= maxiter
        # Computation
        #qx         .-=   dt./(ρ.*dc + dt).*(qx + dc.*diff(C)./dx)                  # Ludo's formulation
        qx          .=   1.0./(dτ + ρ.*dc) .* (ρ.*dc.*qx - dc.*dτ.*diff(C)./dx) # My formulation 
        #C[2:end-1] .-=   dτ.* diff(qx)./dx
        C[2:end-1] .-=   dτ./(1 + dτ/ξ) .* ((C[2:end-1] .- Ceq)./ξ .+ diff(qx)./dx)
        # Visualisation
        if iter % ncheck == 0
            err = maximum(abs.(diff(dc.*diff(C)./dx)./dx .- (C[2:end-1].-Ceq)/ξ))
            push!(iter_evo, iter/ncx); push!(err_evo, err)
            sleep(0.1)
            empty!(ax1)
            empty!(ax2)
            lines!(ax1, xc, C, color = :blue)
            lines!(ax1, xc, C_ini, color = :orange)
            scatter!(ax2, iter_evo, err_evo, color = :blue)
            lines!(ax2, iter_evo, err_evo, color = :blue)
            display(fig1)
        end
        iter += 1
    end        
end

# Call function
diffusionReaction_implicit_1D()