# Solving 2D porous convection equation with central finite differences. Fluid pressure time integration scheme is implicit, while temperature is treated explicitly. Boundary conditions are defined on fluxes.
using CairoMakie, Printf, FFMPEG
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)
# Define Function
@views function porousConvection_2D()
    # Physics
    lx, ly     = 40.0, 20.0            # Length in x and y
    k_ηf       = 1.0                   # Diffusion coefficient
    αρgx,αρgy  = 0.0,1.0
    αρg        = sqrt(αρgx^2+αρgy^2)
    ΔT         = 200.0
    ϕ          = 0.1
    Ra         = 1000
    λ_ρCp      = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    T_Bot      = ΔT/2.0                # BC at domain bottom
    T_Top      = -ΔT/2.0                #BC at domain top
    # Numerics
    re         = 2π
    CFL        = 1.0/sqrt(2.1)
    ncx, ncy   = 127,63               # Number of cells in x and y
    ϵtol       = 1e-8                  # Residual tolerance
    maxiter    = 20*max(ncx,ncy)       # Maximum no. iterations
    ncheck     = ceil(Int,0.25ncx)     # Convergence check frequency
    nvis       = 5                     # Visualisation frequency
    qstp       = 4                     # Quiver density: higher values = less dense
    aspRat     = lx/ly                 # Model aspect ratio
    printFig   = false                 # Printing switch
    writeMovie = false
    # Derived Numerics
    dx        = lx/ncx                # Spatial step size in x
    dy        = ly/ncy                # Spatial step size in y
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1 # Thermal diffusion time step
    θ_dτ      = max(lx,ly)/re/CFL/min(dx,dy)
    β_dτ      = (re*k_ηf)/(CFL*min(dx,dy)*max(lx,ly))
    nt        = 1000 #ncx^2/5              # No. of time steps to compute
    xc        = LinRange(-lx/2.0+dx/2, lx/2.0-dx/2, ncx) # Coordinate array
    yc        = LinRange(-ly+dy/2,         0.0, ncy) # Coordinate array
    # Initialisation
    T         = @. ΔT*exp(-xc^2 - (yc'+ly/2.0)^2); T_ini = copy(T);
    Tx_avg    = zeros(Float64, ncx-1, ncy  )
    Ty_avg    = zeros(Float64, ncx  , ncy-1)
    Pf        = zeros(Float64, ncx  , ncy  )
    qDx       = zeros(Float64, ncx+1, ncy  )
    qDy       = zeros(Float64, ncx  , ncy+1)
    qDx_avg   = zeros(Float64, ncx  , ncy  )
    qDy_avg   = zeros(Float64, ncx  , ncy  )
    qTx       = zeros(Float64, ncx+1, ncy  )
    qTy       = zeros(Float64, ncx  , ncy+1)
    r_Pf      = zeros(Float64, ncx  , ncy  )
    # Plot initial configuration
    fig1      = Figure(figure_padding=30)                 # Plotting
    ax1       = Axis(fig1[1, 1], title="Porous convection 2D", aspect=aspRat, xlabel=L"\textit{xc}", ylabel=L"\textit{yc}", limits=(minimum(xc), maximum(xc), minimum(yc), maximum(yc)))
    heatmap!(ax1, xc, yc, T_ini', colormap=:heat)
    Colorbar(fig1[2,1], vertical=false, colormap=:heat, label=L"\textit{T}", limits=(-ΔT/2, ΔT/2))
    fig1
    # Time loop
    for it=1:nt
        # Calculate fluid pressure (implicit scheme)
        iter = 1; err_Pf = 2ϵtol;
        while err_Pf >= ϵtol && iter <= maxiter
            # Computation of fluid pressure flux
            Tx_avg .= (T[2:end,:] .+ T[1:end-1,:])./2.0; Ty_avg .= (T[:,2:end] .+ T[:,1:end-1])./2.0
            qDx[2:end-1,:].-=   1.0./(1.0 + θ_dτ).*(qDx[2:end-1,:] + k_ηf.*diff(Pf, dims=1)./dx .- αρgx.*Tx_avg) 
            qDy[:,2:end-1].-=   1.0./(1.0 + θ_dτ).*(qDy[:,2:end-1] + k_ηf.*diff(Pf, dims=2)./dy .- αρgy.*Ty_avg)
            # Fluid pressure update
            Pf            .-=   (diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy)./β_dτ
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
        # Boundary conditions
        qTx[[1,end],:] .= 0.0;                              # Left and right
        qTy[:,1  ] .= -λ_ρCp.*(2*(T[:,1  ].-T_Bot))./dy     # Bottom
        qTy[:,end] .= -λ_ρCp.*(2*(T_Top.-T[:,end]))./dy     # Top
        # Compute thermal fluxes
        qTx[2:end-1,:] .= .-λ_ρCp.*diff(T, dims=1)./dx
        qTy[:,2:end-1] .= .-λ_ρCp.*diff(T, dims=2)./dy
        # Temperature update for diffusion
        T  .-= dt.*(diff(qTx, dims=1)./dx + diff(qTy, dims=2)./dy)
        # Temperature update for advection
        T[2:end  ,:] .-= dt./ϕ.*max.(0.0,qDx_avg[2:end,:]  ).*diff(T, dims=1)./dx
        T[1:end-1,:] .-= dt./ϕ.*min.(qDx_avg[1:end-1,:],0.0).*diff(T, dims=1)./dx
        T[:,  2:end] .-= dt./ϕ.*max.(0.0,qDy_avg[:,2:end]  ).*diff(T, dims=2)./dy
        T[:,1:end-1] .-= dt./ϕ.*min.(qDy_avg[:,1:end-1],0.0).*diff(T, dims=2)./dy
        # Visualisation
        if it % nvis == 0
            qDx_avg .= (qDx[2:end,:] + qDx[1:end-1,:])./2.0; qDy_avg .= (qDy[:,2:end] + qDy[:,1:end-1])./2.0
            sleep(1.0)
            empty!(ax1)
            heatmap!(ax1, xc, yc, T, colormap=:heat)
            arrows!(ax1, xc[1:qstp:end], yc[1:qstp:end], qDx_avg[1:qstp:end,1:qstp:end], qDy_avg[1:qstp:end,1:qstp:end], arrowsize=10, lengthscale=0.5, normalize=true)
            display(fig1)
            # Print figure
            if printFig
                save(@sprintf("./ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/doc/png/%04d_BConFlux.png",it),fig1)
            end
        end
    end

    # Write movie
    if writeMovie
        FFMPEG.ffmpeg_exe(`-framerate 12 -f image2 -pattern_type glob -i /Users/lcandiot/Developer/Julia/BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/doc/png/'*'_BConFlux.png -vf "scale=1920:1080" -c:v libx264 -pix_fmt yuv420p -y "/Users/lcandiot/Developer/Julia/BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/doc/out_BConFlux_movie.mov"`)    
    end
end

# Call function
porousConvection_2D()