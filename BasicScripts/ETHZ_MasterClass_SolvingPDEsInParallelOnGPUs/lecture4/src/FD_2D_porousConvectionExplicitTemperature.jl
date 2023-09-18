# Solving 2D porous convection equation with central finite differences
using CairoMakie, Printf
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
# Define Function
@views function porousConvection_2D()
    # Physics
    lx, ly   = 40.0, 20.0                  # Length in x and y
    k_ηf     = 1.0                   # Diffusion coefficient
    # Numerics
    re       = 2π
    CFL      = 1.0/sqrt(2.1)
    ncx, ncy = 400,200              # Number of cells in x and y
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                 # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)     # Convergence check frequency
    aspRat  = lx/ly                 # Model aspect ratio
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size in x
    dy      = ly/ncy                # Spatial step size in y
    θ_dτ    = max(lx,ly)/re/CFL/min(dx,dy)
    β_dτ    = (re*k_ηf)/(CFL*min(dx,dy)*max(lx,ly))
    nt      = 10 #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    yc      = LinRange(dy/2, ly-dy/2, ncy) # Coordinate array
    # Initialisation
    Pf      = @. exp(-(xc-lx/4)^2 -(yc'-ly/4)^2); Pf_ini = copy(Pf);
    qDx     = zeros(Float64, ncx+1, ncy  )
    qDy     = zeros(Float64, ncx  , ncy+1)
    r_Pf    = zeros(Float64, ncx  , ncy  )
    fig1    = Figure()                 # Plotting
    ax1     = Axis(fig1[1, 1], aspect=aspRat, xlabel=L"\textit{xc}", ylabel=L"\textit{yc}")
    heatmap!(ax1, xc, yc, Pf, colormap=:viridis)
    Colorbar(fig1[2,1], vertical=false, colormap=:viridis, label=L"\textit{P}_{\textit{f}}")
    fig1
    # Time loop
    for it=1:nt
        # Iteration loop
        iter = 1; err_Pf = 2ϵtol;
        while err_Pf >= ϵtol && iter <= maxiter
            # Computation
            qDx[2:end-1,:].-=   1.0./(1.0 + θ_dτ).*(qDx[2:end-1,:] + k_ηf.*diff(Pf, dims=1)./dx)                  # Ludo's formulation 
            qDy[:,2:end-1].-=   1.0./(1.0 + θ_dτ).*(qDy[:,2:end-1] + k_ηf.*diff(Pf, dims=2)./dy)                  # Ludo's formulation
            Pf          .-=   (diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy)./β_dτ
            # Convergence check
            if iter % ncheck == 0
                r_Pf .= diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy 
                err_Pf = maximum(abs.(r_Pf))
                # Print convergence to stdout
                @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/ncx,err_Pf)
            end
            iter += 1
        end
        # Print to stdout
        @printf("it = %d, iter/nx=%.1f, err_Pf=%1.3e\n",it,iter/ncx,err_Pf)
        # Visualisation
        sleep(0.1)
        empty!(ax1)
        heatmap!(ax1, xc, yc, Pf, colormap=:viridis)
        Colorbar(fig1[2,1], vertical=false, colormap=:viridis, label=L"\textit{P}_{\textit{f}}")
        fig1
    end
end

# Call function
porousConvection_2D()