using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools, HDF5
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
# Macros
macro d_xa(A) esc(:( $A[ix+1,iy  ] - $A[ix,iy] ))end
macro d_ya(A) esc(:( $A[ix  ,iy+1] - $A[ix,iy] ))end

function Pf_diffusion_2D(; do_check=true, writeOut=true)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    nx = ny     = 16 * 2 .^(1:8)               # Array resolution#nx,ny   = 511,511
    T_eff       = zeros(Float64, length(nx))
    for iRes in eachindex(nx)
        ϵtol    = 1e-8
        maxiter = 500;#max(nx,ny)
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
            # Compute fluxes
            @inbounds compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
            # Update pressure
            @inbounds update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
            if do_check && iter%ncheck == 0
                r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
                err_Pf = maximum(abs.(r_Pf))
                @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx[iRes],err_Pf)
                display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
            end
            iter += 1
        end
        # Monitor timing
        #t_toc = Base.time()                         # End time [s]
        t_toc = @belapsed begin compute!($Pf, $qDx, $qDy, $_dx, $_dy, $_β_dτ, $k_ηf_dx, $k_ηf_dy, $_1_θ_dτ) end
        niter = iter
        A_eff = (3*8)*(nx[iRes]-1)*(ny[iRes]) + (3*8)*(nx[iRes])*(ny[iRes]-1) + (3*8)*(nx[iRes])*(ny[iRes]) # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
        #t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
        t_it = t_toc
        T_eff[iRes] = A_eff/t_it                          # Effective memory throughput [GB/s]
        @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff[iRes]/1e9)
        #@printf("Elapsed time = %1.3e s\n", t_toc)
    end
    # Write output
    if writeOut
        h5open("./data/Pf_diffusion_2D_perf_loop_fun.h5", "w") do outFile
            write(outFile, "monitor/T_eff", T_eff)
            write(outFile, "model/nx",      nx   )
            write(outFile, "model/ny",      ny   ) 
        end
    end
    # Return
    return
end

# Compute functions
function compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
    nx,ny = size(Pf)
    # Compute qDx
    Threads.@threads for iy=1:ny
        for ix=1:nx-1
            qDx[ix+1,iy] -= (qDx[ix+1,iy] + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
        end
    end
    # Compute qDy
    Threads.@threads for iy=1:ny-1
        for ix=1:nx
            qDy[ix,iy+1] -= (qDy[ix,iy+1] + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
        end
    end
    return nothing
end

function update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
    nx,ny = size(Pf)
    Threads.@threads for iy=1:ny
        for ix=1:nx
            Pf[ix,iy] -= (@d_xa(qDx)*_dx + @d_ya(qDy)*_dy)*_β_dτ
        end
    end
    return nothing
end

function compute!(Pf, qDx, qDy, _dx, _dy, _β_dτ, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
   @inbounds compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
   @inbounds update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
   return nothing 
end

Pf_diffusion_2D(; do_check=false, writeOut=false)
