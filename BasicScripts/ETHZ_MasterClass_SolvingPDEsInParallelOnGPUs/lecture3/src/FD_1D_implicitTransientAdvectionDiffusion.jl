# Solving 1D diffusion advection equation with central finite differences and implicit time integration
using Pkg, CairoMakie
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
# Define Function
@views function advectionDiffusion_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    time    = 0.0                   # Time
    vx      = 1.0                   # Fixed free stream velocity
    # Numerics
    ncx     = 200                   # Number of cells in x
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                  # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)      # Convergence check frequency
    nvis    = 1                     # Plotting frequency
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx/abs(vx)            # Time step
    da      = lx^2/dc/dt            # Damköhler number
    re      = π + sqrt(π^2 + da)    # Reynolds number
    ρ       = (lx/(dc*re))^2        # Density
    dτ      = dx/sqrt(1/ρ) #dx^2/dc/2.1           # Pseudo-time step size
    nt      = 10              # No. of time steps to compute
    # Initialisation
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C); C_old = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax1     = Axis(fig1[1, 1], xlabel=L"\textit{x}[m]", ylabel=L"\textit{C}[mol]", limits=(0, lx, 0, 2), title="1D Steady-state Diffusion (impl.)")
    ax2     = Axis(fig1[2, 1], yscale = log10, xlabel=L"\textit{iter/ncx}", ylabel=L"||\textit{r}||_{∞}")
    lines!(ax1, xc, C, color = :blue)
    lines!(ax1, xc, C_ini, color = :orange)
    display(fig1) 
    # TIME LOOP
    for it=1:nt
        # Iteration loop
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err >= ϵtol && iter <= maxiter
            # Compute diffusion
            qx         .-=   dτ./(ρ.*dc + dτ).*(qx + dc.*diff(C)./dx)                  # Ludo's formulation
            # qx          .=   1.0./(dτ + ρ.*dc) .* (ρ.*dc.*qx - dc.*dτ.*diff(C)./dx) # My formulation 
            # C[2:end-1] .-=   dτ.* diff(qx)./dx
            C[2:end-1] .-=   dτ./(1 + dτ/dt) .* ((C[2:end-1] .- C_old[2:end-1])./dt .+ diff(qx)./dx)
            # Convergence check
            if iter % ncheck == 0
                err = maximum(abs.(diff(dc.*diff(C)./dx)./dx .- (C[2:end-1].-C_old[2:end-1])/dt))
                push!(iter_evo, iter/ncx); push!(err_evo, err)
            end
            iter += 1
        end
        # Compute advection
        C[2:end  ] .-= dt.*max(0.0,vx).*diff(C)./dx
        C[1:end-1] .-= dt.*min(vx,0.0).*diff(C)./dx 
        # Visualisation
        if it % nvis == 0
            sleep(0.05)
            empty!(ax1)
            empty!(ax2, )
            lines!(ax1, xc, C, color = :blue)
            ax1.title = L"\textit{t} = %$(time)"
            lines!(ax1, xc, C_ini, color = :orange)
            lines!(ax2, iter_evo, err_evo, color = :blue)
            display(fig1)
        end
        # Update physics
        C_old .= C
        time  += dt
    end        
end

# Call function
advectionDiffusion_implicit_1D()