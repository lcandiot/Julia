# Solving 1D nonlinear diffusion equation with central finite differences
using GLMakie
GLMakie.activate!()
# Define Function
@views function diffusionNonlinear_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    n       = 4                     # Power-law exponent
    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 100                     # Visualise every XXX
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx^2/dc/10           # Time step size
    nt      = ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. 0.5cos(9π*xc/lx)+0.5; C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1])
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Computation
        qx          .= .-dc.*diff(C.^n )./dx
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
diffusionNonlinear_1D()