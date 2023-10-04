# Solving 1D diffusion equation with central finite differences
using Pkg, CairoMakie
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
# Define Function
@views function diffusionAdvection_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    vx      = 10.0                   # Constant velocity
    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 2                     # Visualise every XXX
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt_dif  = dx^2/dc/2.1           # Time step size for diffusion
    dt_adv  = dx/abs(vx)            # Time step size for advection
    dt      = min(dt_dif, dt_adv)   # Minimum time step 
    nt      = ncx^2/100             # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. exp(-(xc-lx/4)^2); C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1], xlabel=L"\textit{x}[m]", ylabel=L"\textit{C}[mol]", limits=(0, lx, 0, 1), title="1D Diffusion-Advection")
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Computation
        qx          .= .-dc.*diff(C )./dx
        C[2:end-1] .-=   dt.*diff(qx)./dx
        # Advection 
        vx>0 ? C[2:end] .-= dt.*vx.*diff(C)/dx : C[1:end-1] .-= dt.vx.diff(C)/dx 
        # Visualisation
        if it%nVis == 0
            sleep(0.05)
            empty!(ax)
            lines!(ax, xc, C, color=:blue)
            lines!(ax, xc, C_ini, color=:orange)
            display(fig1)
        end
    end        
end

# Call function
diffusionAdvection_1D()