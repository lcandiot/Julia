# Solving 1D damped wave equation with central finite differences
using GLMakie
GLMakie.activate!()
# Define Function
@views function dampedWaveEquation_explicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    ρ       = 20.0                  # Density
    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 50                     # Visualise every XXX
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx^2/dc/2.1           # Time step size
    nt      = ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1])
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Computation
        qx         .-=   dt.*(diff(C)/dx + qx/dc)/ρ
        C[2:end-1] .-=   dt.*diff(qx)./dx
        # Visualisation
        if it%nVis == 0
            sleep(0.1)
            fig1    = Figure()                 # Plotting
            ax      = Axis(fig1[1, 1])
            lines!(ax, xc, C)
            lines!(ax, xc, C_ini)
            display(fig1)
            DataInspector(fig1)
        end
    end        
end

# Call function
dampedWaveEquation_explicit_1D()