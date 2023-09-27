using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
# Macros
macro d_xa(A) esc(:( $A[ix+1,iy  ] - $A[ix,iy] ))end
macro d_ya(A) esc(:( $A[ix  ,iy+1] - $A[ix,iy] ))end

function Pf_diffusion_2D(; do_check=true)
    # physics
    lx,ly   = 20.0,20.0
    k_ηf    = 1.0
    # numerics
    ncx,ncy   = 511,511
    ϵtol    = 1e-8
    maxiter = max(ncx,ncy)
    ncheck  = ceil(Int,0.25max(ncx,ncy))
    cfl     = 1.0/sqrt(2.1)
    re      = 2π
    # derived numerics
    dx,dy   = lx/ncx,ly/ncy
    xc,yc   = LinRange(dx/2,lx-dx/2,ncx),LinRange(dy/2,ly-dy/2,ncy)
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
    qDx,qDy = zeros(Float64, ncx+1,ncy),zeros(Float64, ncx,ncy+1)
    r_Pf    = zeros(ncx,ncy)
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
            @printf("  iter/ncx=%.1f, err_Pf=%1.3e\n",iter/ncx,err_Pf)
            display(heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1
    end
    # Monitor timing
    #t_toc = Base.time()                         # End time [s]
    t_toc = @belapsed compute!($Pf, $qDx, $qDy, $_dx, $_dy, $_β_dτ, $k_ηf_dx, $k_ηf_dy, $_1_θ_dτ)
    niter = iter
    A_eff = (3*8)*(ncx-1)*(ncy) + (3*8)*(ncx)*(ncy-1) + (3*8)*(ncx)*(ncy) # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
    t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
    T_eff = A_eff/t_it                          # Effective memory throughput [GB/s]
    # @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff/1e9)
    @printf("Elapsed time = %1.3e s\n", t_toc)
    # Monitor one iteration
    
    return
end

# Compute functions
function compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
    ncx,ncy = size(Pf)
    # Compute qDx
    for iy=1:ncy
        for ix=1:ncx-1
            qDx[ix+1,iy] -= (qDx[ix+1,iy] + k_ηf_dx*@d_xa(Pf))*_1_θ_dτ
        end
    end
    # Compute qDy
    for iy=1:ncy-1
        for ix=1:ncx
            qDy[ix,iy+1] -= (qDy[ix,iy+1] + k_ηf_dy*@d_ya(Pf))*_1_θ_dτ
        end
    end
    return nothing
end

function update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
    ncx,ncy = size(Pf)
    for iy=1:ncy
        for ix=1:ncx
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

Pf_diffusion_2D(; do_check=false)
