using Plots,Plots.Measures,Printf, CairoMakie, BenchmarkTools, HDF5
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

function Pf_diffusion_2D(; do_check=true, writeOut=true)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    nx = ny     = 16 * 2 .^(1:8)               # Array resolution#nx,ny   = 511,511
    T_eff       = zeros(Float64, length(nx))
    for iRes in eachindex(nx)
        ϵtol    = 1e-8
        maxiter = 500
        ncheck  = ceil(Int,0.25max(nx[iRes],ny[iRes]))
        cfl     = 1.0/sqrt(2.1)
        re      = 2π
        # derived numerics
        dx,dy   = lx/nx[iRes],ly/ny[iRes]
        xc,yc   = LinRange(dx/2,lx-dx/2,nx[iRes]),LinRange(dy/2,ly-dy/2,ny[iRes])
        θ_dτ    = max(lx,ly)/re/cfl/min(dx,dy)
        β_dτ    = (re*k_ηf)/(cfl*min(dx,dy)*max(lx,ly))
        # inverse multiplications
        k_ηf_dx, k_ηf_dy = k_ηf/dx, k_ηf/dy
        _1_θ_dτ = 1.0/(1.0 + θ_dτ)
        _dx     = 1.0/dx
        _dy     = 1.0/dy
        _β_dτ   = 1.0/β_dτ
        # array initialisation
        Pf      = @. exp(-(xc-lx/2)^2 -(yc'-ly/2)^2)
        qDx,qDy = zeros(Float64, nx[iRes]+1,ny[iRes]),zeros(Float64, nx[iRes],ny[iRes]+1)
        r_Pf    = zeros(nx[iRes],ny[iRes])
        # iteration loop
        iter = 1; err_Pf = 2ϵtol; t_tic = 0; warmup_iter = 11
        while err_Pf >= ϵtol && iter <= maxiter
            # Define starting time
            if iter==warmup_iter
                t_tic = Base.time()
            end
            qDx[2:end-1,:] .-= (qDx[2:end-1,:]        .+ (diff(Pf,dims=1).*k_ηf_dx)).*_1_θ_dτ
            qDy[:,2:end-1] .-= (qDy[:,2:end-1]        .+ (diff(Pf,dims=2).*k_ηf_dy)).*_1_θ_dτ
            Pf             .-= (diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy     ).*_β_dτ
            #compute_flux!(qDx, qDy, Pf, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
            #update_pressure!(qDx, qDy, Pf, _dx, _dy, _β_dτ)
            if do_check && iter%ncheck == 0
                r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
                err_Pf = maximum(abs.(r_Pf))
                @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx[iRes],err_Pf)
                display(Plots.heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
            end
            iter += 1
        end
        # Monitor timing
        t_toc = Base.time()                         # End time [s]
        niter = iter
        A_eff = (3*8)*(nx[iRes]-1)*(ny[iRes]) + (3*8)*(nx[iRes])*(ny[iRes]-1) + (3*8)*(nx[iRes])*(ny[iRes])          # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
        t_it  =  @belapsed compute!($qDx, $qDy, $Pf, $k_ηf_dx, $k_ηf_dy, $_1_θ_dτ, $_dx, $_dy, $_β_dτ)
        #t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
        T_eff[iRes] = A_eff/t_it                          # Effective memory throughput [GB/s]
        @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff[iRes]/1e9)
    end
    # Write output
    if writeOut
        h5open("./data/Pf_diffusion_2D_perf.h5", "w") do outFile
            write(outFile, "monitor/T_eff", T_eff)
            write(outFile, "model/nx",      nx   )
            write(outFile, "model/ny",      ny   ) 
        end
    end
    # Return
    return
end

function compute_flux!(qDx, qDy, Pf, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
    qDx[2:end-1,:] .-= (qDx[2:end-1,:]        .+ (diff(Pf,dims=1).*k_ηf_dx)).*_1_θ_dτ
    qDy[:,2:end-1] .-= (qDy[:,2:end-1]        .+ (diff(Pf,dims=2).*k_ηf_dy)).*_1_θ_dτ
    # Return
    return nothing
end

function update_pressure!(qDx, qDy, Pf, _dx, _dy, _β_dτ)
    Pf             .-= (diff(qDx,dims=1).*_dx .+ diff(qDy,dims=2).*_dy     ).*_β_dτ
    # Return
    return nothing
end

function compute!(qDx, qDy, Pf, k_ηf_dx, k_ηf_dy, _1_θ_dτ, _dx, _dy, _β_dτ)
    @inbounds compute_flux!(qDx, qDy, Pf, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
    @inbounds update_pressure!(qDx, qDy, Pf, _dx, _dy, _β_dτ)
    # Return
    return nothing
end

Pf_diffusion_2D(; do_check=false, writeOut=true)
