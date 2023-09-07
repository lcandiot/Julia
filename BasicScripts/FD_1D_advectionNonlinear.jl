# Solving 1D nonlinear Advection equation with central finite differences
using GLMakie
GLMakie.activate!()
# Define Function
@views function advectionNonlinear_1D()
    # Physics
    lx        = 20.0                  # Length in x
    vx        = 1.0                  # Constant velocity
    n         = 2                     # Power-law exponent
    ttot      = 20.0                    # Total time
    ttot_half = ttot/2.0            # Half time of simulation
    time      = 0.0                   # Current time
    # Numerics
    ncx       = 1000                   # Number of cells in x
    nVis      = 15                     # Visualise every XXX
    vel_switch = 1                  
    # Derived Numerics
    dx        = lx/ncx                # Spatial step size
    dt        = dx/abs(vx)/2            # Time step size for advection
    nt        = 2ncx             # No. of time steps to compute
    xc        = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    C         = @. exp(-(xc-lx/4)^2); C_ini = copy(C)
    fig1      = Figure()                 # Plotting
    ax        = Axis(fig1[1, 1])
    lines!(ax, xc, C)
    lines!(ax, xc, C_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Update time
        time = time .+ dt
        # Advection 
        C[2:end  ] .-= dt.*max(vx,0.0).*diff(C.^n)./dx
        C[1:end-1] .-= dt.*min(vx,0.0).*diff(C.^n)./dx
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
        # Switch velocity at half time
        if time >= ttot_half && vel_switch == 1
            vx = -vx
            vel_switch = 0
        end
    end        
end

# Call function
advectionNonlinear_1D()