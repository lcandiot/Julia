# Solving 2D porous convection equation with central finite differences
using CairoMakie, Printf
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
# Define Function
@views function porousConvection_2D()
    # Physics
    lx, ly    = 40.0, 20.0            # Length in x and y
    k_ηf      = 1.0                   # Diffusion coefficient
    αρgx,αρgy = 0.0,1.0
    αρg       = sqrt(αρgx^2+αρgy^2)
    ΔT        = 200.0
    ϕ         = 0.1
    Ra        = 100
    λ_ρCp     = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    # Numerics
    re        = 2π
    CFL       = 1.0/sqrt(2.1)
    ncx, ncy  = 400,200               # Number of cells in x and y
    ϵtol      = 1e-8                  # Residual tolerance
    maxiter   = 20*max(ncx,ncy)       # Maximum no. iterations
    ncheck    = ceil(Int,0.25ncx)     # Convergence check frequency
    aspRat    = lx/ly                 # Model aspect ratio
    # Derived Numerics
    dx        = lx/ncx                # Spatial step size in x
    dy        = ly/ncy                # Spatial step size in y
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1 # Thermal diffusion time step
    θ_dτ      = max(lx,ly)/re/CFL/min(dx,dy)
    β_dτ      = (re*k_ηf)/(CFL*min(dx,dy)*max(lx,ly))
    nt        = 10 #ncx^2/5              # No. of time steps to compute
    xc        = LinRange(-lx/2.0+dx/2, lx/2.0-dx/2, ncx) # Coordinate array
    yc        = LinRange(-ly/2.0+dy/2, ly/2.0-dy/2, ncy) # Coordinate array
    # Initialisation
    T         = @. ΔT*exp(-xc^2 - (yc')^2); T_ini = copy(T);
    Pf        = zeros(Float64, ncx  , ncy)
    qDx       = zeros(Float64, ncx+1, ncy  )
    qDy       = zeros(Float64, ncx  , ncy+1)
    qTx       = zeros(Float64, ncx-1, ncy  )
    qTy       = zeros(Float64, ncx  , ncy-1)
    r_Pf      = zeros(Float64, ncx  , ncy  )
    fig1      = Figure()                 # Plotting
    ax1       = Axis(fig1[1, 1], aspect=aspRat, xlabel=L"\textit{xc}", ylabel=L"\textit{yc}")
    ax2       = Axis(fig1[2,1])
    heatmap!(ax1, xc, yc, T_ini, colormap=:heat)
    #Colorbar(fig, vertical=false, colormap=:heat, label=L"\textit{T}")
    fig1
    # Time loop
    for it=1:nt
        # Calculate fluid pressure (implicit scheme)
        iter = 1; err_Pf = 2ϵtol;
        while err_Pf >= ϵtol && iter <= maxiter
            # Computation of fluid pressure flux
            qDx[2:end-1,:].-=   1.0./(1.0 + θ_dτ).*(qDx[2:end-1,:] + k_ηf.*diff(Pf, dims=1)./dx)                  # Ludo's formulation 
            qDy[:,2:end-1].-=   1.0./(1.0 + θ_dτ).*(qDy[:,2:end-1] + k_ηf.*diff(Pf, dims=2)./dy)                  # Ludo's formulation
            # Fluid pressure update
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
        # Calculate temperature (explicit scheme)
        dt_adv = ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1
        dt     = min(dt_diff,dt_adv) # Time step calculation
        # Compute thermal fluxes
        qTx .= .-λ_ρCp.*diff(T, dims=1)./dx
        qTy .= .-λ_ρCp.*diff(T, dims=2)./dy
        # Temperture update for diffusion
        T[2:end-1, 2:end-1] .-= dt.*(diff(qTx[:,2:end-1], dims=1)./dx + diff(qTy[2:end-1,:], dims=2)./dy)
        # Temperature update for advection
        T[2:end  ,:] .-= dt./ϕ.*max.(0.0,qDx[2:end-1,:]).*diff(T, dims=1)./dx
        T[1:end-1,:] .-= dt./ϕ.*min.(qDx[1:end-2,:],0.0).*diff(T, dims=1)./dx
        T[:,  2:end] .-= dt./ϕ.*max.(0.0,qDy[:,2:end-1]).*diff(T, dims=2)./dy
        T[:,1:end-1] .-= dt./ϕ.*min.(qDy[:,1:end-2],0.0).*diff(T, dims=2)./dy
        # Boundary conditions
        T[:,1] .= ΔT/2; T[:,end] .= -ΔT/2   # Top and bottom
        T[[1,end],:] .= T[[2,end-1],:]      # Left and right
        # Visualisation
        sleep(1.0)
        empty!(ax1)
        heatmap!(ax1, xc, yc, T, colormap=:heat)
        Colorbar(fig1[2,1], vertical=false, colormap=:heat, label=L"\textit{T}", limits=(-ΔT/2, ΔT/2))
        display(fig1)
    end
end

# Call function
porousConvection_2D()