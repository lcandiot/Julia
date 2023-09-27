using Plots,Plots.Measures,Printf
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

function Pf_diffusion_2D(; do_check=true)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    nx,ny   = 511,511
    ϵtol    = 1e-8
    maxiter = max(nx,ny)
    ncheck  = ceil(Int,0.25max(nx,ny))
    cfl     = 1.0/sqrt(2.1)
    re      = 2π
    # derived numerics
    dx,dy   = lx/nx,ly/ny
    xc,yc   = LinRange(dx/2,lx-dx/2,nx),LinRange(dy/2,ly-dy/2,ny)
    θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
    β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
    # array initialisation
    Pf      = @. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2)
    qDx,qDy = zeros(Float64, nx+1,ny),zeros(Float64, nx,ny+1)
    r_Pf    = zeros(nx,ny)
    # iteration loop
    iter = 1; err_Pf = 2ϵtol; t_tic = 0; warmup_iter = 11
    while err_Pf >= ϵtol && iter <= maxiter
        # Define starting time
        if iter==warmup_iter
            t_tic = Base.time()
        end
        qDx[2:end-1,:] .-= (qDx[2:end-1,:] .+ k_ηf.*(diff(Pf,dims=1)./dx))./(1.0 + θ_dτ)
        qDy[:,2:end-1] .-= (qDy[:,2:end-1] .+ k_ηf.*(diff(Pf,dims=2)./dy))./(1.0 + θ_dτ)
        Pf             .-= (diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy)./β_dτ
        if do_check && iter%ncheck == 0
            r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1
    end
    # Monitor timing
    t_toc = Base.time()                         # End time [s]
    niter = iter
    A_eff = (3*8)*(nx-1)*(ny) + (3*8)*(nx)*(ny-1) + (3*8)*(nx)*(ny)          # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
    t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
    T_eff = A_eff/t_it                          # Effective memory throughput [GB/s]
    @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff/1e9)
    return
end

Pf_diffusion_2D(; do_check=false)
