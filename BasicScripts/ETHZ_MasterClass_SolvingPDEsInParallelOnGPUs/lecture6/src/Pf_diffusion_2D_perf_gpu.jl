# Running 2D fluid pressure diffusion on Apple Silicon GPU
using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools, HDF5, Metal
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

# Macros
macro d_xa(A) esc(:( $A[ix+1,iy  ] - $A[ix,iy] ))end
macro d_ya(A) esc(:( $A[ix  ,iy+1] - $A[ix,iy] ))end

function Pf_diffusion_2D(; do_check=true, writeOut=true, writeTest=true)
    # Physics
    lx,ly   = Float32(20.0), Float32(20.0)
    k_ηf    = Float32(1.0)
    
    # Numerics
    nx = ny = Int32(127)               # Array resolution#nx,ny   = 511,511
    ϵtol    = Float32(1e-7)
    maxiter = Int32(1e4);#max(nx,ny)
    ncheck  = Int32( ceil(Int32, 0.25max(nx, ny)) )
    cfl     = Float32( 1.0 / sqrt(2.1) )
    re      = Float32(2π)
    
    # GPU device
    MTLDevice(1)
    threads = (16,16) # 16 by 16 seems to be the most performant configuration
    groups  = cld.((nx, ny), threads)

    # Derived numerics
    dx,dy   = Float32( lx / nx ),Float32( ly / ny )
    xc,yc   = LinRange{Float32}(dx/2, lx - dx/2, nx), LinRange{Float32}(dy/2, ly - dy/2, ny)
    θ_dτ    = Float32( max(lx,ly) / re / cfl / min(dx,dy))
    β_dτ    = Float32( (re * k_ηf) / (cfl * min(dx,dy) * max(lx,ly)) )
    
    # Inverse multiplications
    k_ηf_dx, k_ηf_dy = Float32(k_ηf/dx), Float32(k_ηf/dy)
    _1_θ_dτ = Float32( 1.0 / (1.0 + θ_dτ) )
    _dx     = Float32( 1.0 / dx )
    _dy     = Float32( 1.0 / dy )
    _β_dτ   = Float32( 1.0 / β_dτ )
    
    # Array initialisation
    Pf   = MtlArray( zeros(Float32, nx  , ny) )
    Pf  .= MtlArray( [exp( -( (ix-1) * dx - lx/Float32(2.0) )^2 - ( (iy-1)*dy - ly/Float32(2.0) )^2 ) 
                        for ix=1:size(Pf,1), iy=1:size(Pf,2)] )
    qDx  = MtlArray( zeros(Float32, nx+1, ny  ) )
    qDy  = MtlArray( zeros(Float32, nx  , ny+1) )
    r_Pf = MtlArray( zeros(Float32, nx  , ny  ) )
    
    # Iteration loop
    iter = Int32(1); err_Pf = Float32(2ϵtol); t_tic = Float32(0.0); warmup_iter = Int32(11)
    while err_Pf >= ϵtol && iter <= maxiter
        # Define starting time
        if ( iter==warmup_iter )
            t_tic = Float32( Base.time() )
        end
        
        # Solve
        compute!(Pf, qDx, qDy, _dx, _dy, _β_dτ, k_ηf_dx, k_ηf_dy, _1_θ_dτ, threads, groups)

        # Check for convergence
        if do_check && iter%ncheck == 0
            r_Pf  .= diff( qDx, dims=1 ) ./ dx .+ diff( qDy, dims=2 ) ./ dy
            err_Pf = maximum( abs.(r_Pf) )
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            Plots.display(Plots.heatmap(xc,yc,Array(Pf)';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1

        # Write pressure field for unit testing
        if (iter==50 && writeTest)
            h5open("./test/gpu_testData.h5", "w") do outfile
                write(outfile, "centers/Pf", Array{Float64}(Pf))
            end
        end
    end
    
    # Monitor timing
    # t_toc = @belapsed begin compute!($Pf, $qDx, $qDy, $_dx, $_dy, $_β_dτ, $k_ηf_dx, $k_ηf_dy, $_1_θ_dτ, $threads, $groups) end
    # niter = iter
    # A_eff = (3*8)*(nx-1)*(ny) + (3*8)*(nx)*(ny-1) + (3*8)*(nx)*(ny) # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
    #t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
    # t_it = t_toc
    # T_eff = A_eff/t_it                          # Effective memory throughput [GB/s]
    # @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff/1e9)
    
    # Return
    return
end

# Compute fluxes
function compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
    nx,ny  = size(Pf)
    ix, iy = thread_position_in_grid_2d()

    # qDx
    if ( ix>=1 && ix<=nx-1 && iy>=1 && iy<=ny )
        @inbounds qDx[ix+1,iy] -= ( qDx[ix+1,iy] + k_ηf_dx * @d_xa(Pf) ) * _1_θ_dτ
    end

    # qDy
    if ( ix>=1 && ix<=nx && iy>=1 && iy<=ny-1 )
        @inbounds qDy[ix,iy+1] -= ( qDy[ix,iy+1] + k_ηf_dy * @d_ya(Pf) ) * _1_θ_dτ
    end

    # Return
    return nothing
end

# Update pressure
function update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
    nx,ny = size(Pf)
    ix, iy = thread_position_in_grid_2d()

    if ( ix>=1 && ix<=nx && iy>=1 && iy<=ny )
        @inbounds Pf[ix,iy] -= ( @d_xa(qDx) * _dx + @d_ya(qDy) * _dy ) * _β_dτ
    end

    # Return
    return nothing
end

# Solver function
function compute!(Pf, qDx, qDy, _dx, _dy, _β_dτ, k_ηf_dx, k_ηf_dy, _1_θ_dτ, threads, groups)
   Metal.@sync @metal threads=threads groups=groups compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
   Metal.@sync @metal threads=threads groups=groups update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
   return nothing 
end

# Run
Pf_diffusion_2D(; do_check=true, writeOut=false, writeTest=false)
