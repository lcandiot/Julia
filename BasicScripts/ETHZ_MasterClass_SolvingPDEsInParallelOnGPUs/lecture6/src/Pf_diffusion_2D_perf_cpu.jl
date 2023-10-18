using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools, HDF5
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
# Macros
macro d_xa(A) esc(:( $A[ix+1,iy  ] - $A[ix,iy] ))end
macro d_ya(A) esc(:( $A[ix  ,iy+1] - $A[ix,iy] ))end

function Pf_diffusion_2D(; do_check=true, writeOut=true, writeTest=true)
    # physics
    lx,ly   = Float32(20.0), Float32(20.0)
    k_ηf    = Float32(1.0)
    # numerics
    nx, ny  = Int32(127), Int32(127)               # Array resolution#nx,ny   = 511,511
    ϵtol    = Float32(1e-7)
    maxiter = Int32(1e4) #max(nx,ny)
    ncheck  = Int32( ceil(Int, 0.25max(nx, ny)) )
    cfl     = Float32(1.0 / sqrt(2.1))
    re      = Float32(2π)
    # derived numerics
    dx,dy   = Float32(lx / nx), Float32(ly / ny)
    xc,yc   = LinRange{Float32}(dx/2, lx - dx/2, nx), LinRange{Float32}(dy/2, ly - dy/2, ny)
    θ_dτ    = Float32( max(lx, ly) / re / cfl / min(dx, dy) )
    β_dτ    = Float32( (re * k_ηf) / (cfl * min(dx,dy) * max(lx,ly)) )
    # inverse multiplications
    k_ηf_dx, k_ηf_dy = Float32( k_ηf / dx ), Float32( k_ηf / dy )
    _1_θ_dτ = Float32( 1.0 / (1.0 + θ_dτ) )
    _dx     = Float32( 1.0 / dx )
    _dy     = Float32( 1.0 / dy )
    _β_dτ   = Float32( 1.0 / β_dτ )
    # array initialisation
    Pf  .= [exp( -((ix-1)*dx - lx/Float32(2.0))^2 - ((iy-1)*dy - ly/Float32(2.0))^2 ) for ix=1:size(Pf,1), iy=1:size(Pf,2)]
    qDx,qDy = zeros(Float32, nx+1, ny), zeros(Float32, nx, ny+1)
    r_Pf    = zeros(Float32, nx  , ny)
    # iteration loop
    iter = Int32(1); err_Pf = Float32(2ϵtol); t_tic = Float32(0.0); warmup_iter = Int32(11)
    while err_Pf >= ϵtol && iter <= maxiter
        # Define starting time
        if iter==warmup_iter
            t_tic = Float32( Base.time() )
        end
        # Compute fluxes
        @inbounds compute_flux!(Pf, qDx, qDy, k_ηf_dx, k_ηf_dy, _1_θ_dτ)
        # Update pressure
        @inbounds update_pressure!(Pf, qDx, qDy, _dx, _dy, _β_dτ)
        if do_check && iter%ncheck == 0
            r_Pf  .= diff(qDx,dims=1)./dx .+ diff(qDy,dims=2)./dy
            err_Pf = maximum(abs.(r_Pf))
            @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx,err_Pf)
            Plots.display(Plots.heatmap(xc,yc,Pf';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        end
        iter += 1

        # Write pressure field for unit testing
        if (iter==50 && writeTest)
            h5open("./test/cpu_testData.h5", "w") do outfile
                write(outfile, "centers/Pf", Pf)
            end
        end
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

Pf_diffusion_2D(; do_check=true, writeOut=false, writeTest=false)
