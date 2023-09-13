# Solving 2D advection diffusion equation with central finite differences
using CairoMakie#GLMakie
CairoMakie.activate!()#GLMakie.activate!()
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
# Define Function
@views function advectionDiffusion_implicit_1D()
    # Physics
    lx, ly  = 10.0, 10.0            # Length in x and y
    dc      = 1.0                   # Diffusion coefficient
    time    = 0.0                   # Time
    vx      = 10.0                  # Free stream velocity x-component
    vy      = -10.0                 # Free stream velocity y-component
    # Numerics
    ncx, ncy = 200, 201              # Number of cells in x and y
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                  # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)      # Convergence check frequency
    nvis    = 1                     # Plotting frequency
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size in x
    dy      = ly/ncy                # Spatial step size in y
    da      = 1000.0                # Damköhler number
    re      = π + sqrt(π^2 + da)    # Reynolds number
    ρ       = (lx/(dc*re))^2        # Density
    dt      = min(dx/abs(vx), dy/abs(vy))/2   # Time step size
    dτ      = min(dx, dy)/sqrt(1/ρ)/sqrt(2.0) # Pseudo-time step size
    nt      = 50              # No. of time steps to compute
    # Initialisation
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array x
    yc      = LinRange(dy/2, ly-dy/2, ncy) # Coordinate array y
    C       = @. exp(-(xc-lx/4)^2 -(yc'-3ly/4)^2); C_old = copy(C)
    qx      = zeros(Float64, ncx-1, ncy  )
    qy      = zeros(Float64, ncx  , ncy-1)
    fig1    = Figure(resolution = (1500, 1500))                 # Plotting
    ax1     = Axis(fig1[1, 1], xlabel=L"\textit{lx}", ylabel=L"\textit{ly}")
    ax2     = Axis(fig1[2, 1:2], yscale = log10, xlabel=L"\textit{iter/nx}", ylabel=L"\textit{err}")
    update_theme!(fontsize=30)
    heatmap!(ax1, xc, yc, C)
    display(fig1) 
    # TIME LOOP
    for it=1:nt
        # Iteration loop
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        while err >= ϵtol && iter <= maxiter
            # Compute diffusion
            qx         .-=   dτ./(ρ.*dc + dτ).*(qx + dc.*diff(C, dims=1)./dx)                  # Ludo's formulation
            qy         .-=   dτ./(ρ.*dc + dτ).*(qy + dc.*diff(C, dims=2)./dy)                  # Ludo's formulation
            C[2:end-1, 2:end-1] .-=   dτ./(1 + dτ/dt) .* ((C[2:end-1, 2:end-1] .- C_old[2:end-1, 2:end-1])./dt .+ diff(qx[:,2:end-1], dims=1)./dx .+ diff(qy[2:end-1,:], dims=2)./dy)
            # Convergence check
            if iter % ncheck == 0
                err = maximum(abs.(diff(dc.*diff(C[:,2:end-1], dims=1)./dx, dims=1)./dx .- (C[2:end-1,2:end-1].-C_old[2:end-1,2:end-1])./dt .+ diff(dc.*diff(C[2:end-1,:], dims=2)./dy, dims=2)./dy))
                push!(iter_evo, iter/ncx); push!(err_evo, err)
            end
            iter += 1
        end
        # Compute advection
        C[2:end  ,:] .-= dt.*max(0.0,vx).*diff(C, dims=1)./dx
        C[1:end-1,:] .-= dt.*min(vx,0.0).*diff(C, dims=1)./dx
        C[:,2:end  ] .-= dt.*max(0.0,vy).*diff(C, dims=2)./dy
        C[:,1:end-1] .-= dt.*min(vy,0.0).*diff(C, dims=2)./dy
        # Visualisation
        if it % nvis == 0
            sleep(0.1)
            empty!(ax1)
            empty!(ax2, )
            heatmap!(ax1, xc, yc, C)
            ax1.title = L"\textit{t} = %$(time)"
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