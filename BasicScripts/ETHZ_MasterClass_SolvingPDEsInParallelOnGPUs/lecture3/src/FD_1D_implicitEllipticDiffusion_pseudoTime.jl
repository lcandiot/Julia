# Solving 1D diffusion equation to steady state with central finite differences
using Pkg, CairoMakie
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
# Define Function
@views function steadyDiffusion_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    ρ       = (lx/(dc*2π))^2        # Numerical density
    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 50                     # Visualise every XXX
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx/sqrt(1/ρ) #dx^2/dc/2.1           # Time step size
    nt      = 5ncx #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1], xlabel=L"\textit{x}[m]", ylabel=L"\textit{C}[mol]", limits=(0, lx, 0, 2), title="1D Steady-state Diffusion (impl.)")
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Computation
        #qx         .-=   dt./(ρ.*dc + dt).*(qx + dc.*diff(C)./dx)                  # Ludo's formulation
        qx          .=   1.0./(dt + ρ.*dc) .* (ρ.*dc.*qx - dc.*dt.*diff(C)./dx) # My formulation 
        C[2:end-1] .-=   dt.* diff(qx)./dx
        # Visualisation
        if it%nVis == 0
            sleep(0.05)
            empty!(ax)
            lines!(ax, xc, C, color = :blue)
            lines!(ax, xc, C_ini, color = :orange)
            display(fig1)
        end
    end        
end

# Call function
steadyDiffusion_implicit_1D()