# Solving 1D reaction-diffusion equation with central finite differences
using GLMakie
GLMakie.activate!()
# Define Function
@views function reactionDiffusion_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 0.1                   # Diffusion coefficient
    ζ       = 10.0                  # Reaction rate
    Ceq     = 0.4                   # Equilibrium concentration
    ttot    = 20.0                  # Total time 
    time    = 0.0                   # Current time
    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 2                     # Visualise every XXX
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx^2/dc/2.1           # Time step size
    nt      = ncx^2/100              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C       = @. exp(-(xc-lx/4)^2); C_ini = copy(C)
    qx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1])
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Update time
        time = time .+ dt
        # Diffusion
        qx           .= .-dc.*diff(C )./dx
        C[2:end-1]  .-=   dt.*diff(qx)./dx
        # Reaction
        C           .-= dt.*(C.-Ceq)./ζ
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
        # Break if simulation time reached end
        # if time >= ttot
        #     break
        # end
    end        
end

# Call function
reactionDiffusion_1D()