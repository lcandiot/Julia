# Benchmarking 2D fluid pressure diffusion on Apple Silicon GPU
using CairoMakie, Printf, Pkg, BenchmarkTools, HDF5, Metal, LaTeXStrings
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

# Macros
macro d_xa(A) esc(:( $A[ix+1,iy  ] - $A[ix,iy] ))end
macro d_ya(A) esc(:( $A[ix  ,iy+1] - $A[ix,iy] ))end

function Pf_diffusion_2D(; do_check=true, writeOut=true, writeTest=true)
    # Physics
    lx,ly   = Float32(20.0), Float32(20.0)
    k_ηf    = Float32(1.0)
    
    # Numerics
    nx = ny = 32 .* 2 .^(0:11)               # Array resolution#nx,ny   = 511,511
    ϵtol    = Float32(1e-7)
    maxiter = Int32(1e4);#max(nx,ny)
    T_peak  = Float64(167.84e9)
    T_apple = Float64(200.0e9)
    
    # Initialise storage arrays for plotting
    T_eff_all   = Float64[ ]
    Res_all     = Float64[ ]
    T_peak_all  = Float64[ ]
    T_apple_all = Float64[ ]

    # Resolution loop
    for iRes in eachindex(nx)
        # Resolution-dependent numerics 
        ncheck  = Int32( ceil(Int32, 0.25max(nx[iRes], ny[iRes])) )
        cfl     = Float32( 1.0 / sqrt(2.1) )
        re      = Float32(2π)
        maxBuff = 11 * nx[iRes] * ny[iRes] * sizeof(Float32)

        # GPU device
        device  = MTLDevice(1)
        threads = (16,16) # 16 by 16 seems to be the most performant configuration
        groups  = cld.((nx[iRes], ny[iRes]), threads)

        # Sanity check, memory allocation, and initialisation
        if (maxBuff > device.maxBufferLength)
            break
        end

        # Derived numerics
        dx,dy   = Float32( lx / nx[iRes] ),Float32( ly / ny[iRes] )
        xc,yc   = LinRange{Float32}(dx/2, lx - dx/2, nx[iRes]), LinRange{Float32}(dy/2, ly - dy/2, ny[iRes])
        θ_dτ    = Float32( max(lx,ly) / re / cfl / min(dx,dy))
        β_dτ    = Float32( (re * k_ηf) / (cfl * min(dx,dy) * max(lx,ly)) )
        
        # Inverse multiplications
        k_ηf_dx, k_ηf_dy = Float32(k_ηf/dx), Float32(k_ηf/dy)
        _1_θ_dτ = Float32( 1.0 / (1.0 + θ_dτ) )
        _dx     = Float32( 1.0 / dx )
        _dy     = Float32( 1.0 / dy )
        _β_dτ   = Float32( 1.0 / β_dτ )
        
        # Array initialisation
        Pf   = MtlArray( zeros(Float32, nx[iRes]  , ny[iRes]) )
        Pf  .= MtlArray( [exp( -( (ix-1) * dx - lx/Float32(2.0) )^2 - ( (iy-1)*dy - ly/Float32(2.0) )^2 ) 
                            for ix=1:size(Pf,1), iy=1:size(Pf,2)] )
        qDx  = MtlArray( zeros(Float32, nx[iRes]+1, ny[iRes]  ) )
        qDy  = MtlArray( zeros(Float32, nx[iRes]  , ny[iRes]+1) )
        # r_Pf = MtlArray( zeros(Float32, nx[iRes]  , ny[iRes]  ) )
        
        # Iteration loop
        # iter = Int32(1); err_Pf = Float32(2ϵtol); t_tic = Float32(0.0); warmup_iter = Int32(11)
        # while err_Pf >= ϵtol && iter <= maxiter
        #     # Define starting time
        #     if ( iter==warmup_iter )
        #         t_tic = Float32( Base.time() )
        #     end
            
        #     # Solve
        #     compute!(Pf, qDx, qDy, _dx, _dy, _β_dτ, k_ηf_dx, k_ηf_dy, _1_θ_dτ, threads, groups)

        #     # Check for convergence
        #     if do_check && iter%ncheck == 0
        #         r_Pf  .= diff( qDx, dims=1 ) ./ dx .+ diff( qDy, dims=2 ) ./ dy
        #         err_Pf = maximum( abs.(r_Pf) )
        #         @printf("  iter/nx=%.1f, err_Pf=%1.3e\n",iter/nx[iRes],err_Pf)
        #         Plots.display(Plots.heatmap(xc,yc,Array(Pf)';xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo,clim=(0,1)))
        #     end
        #     iter += 1

        #     # Write pressure field for unit testing
        #     if (iter==50 && writeTest)
        #         h5open("./test/gpu_testData.h5", "w") do outfile
        #             write(outfile, "centers/Pf", Array{Float64}(Pf))
        #         end
        #     end
        # end
        
        # Monitor timing
        t_it = @belapsed begin compute!($Pf, $qDx, $qDy, $_dx, $_dy, $_β_dτ, $k_ηf_dx, $k_ηf_dy, $_1_θ_dτ, $threads, $groups) end
        # niter = (iter-warmup_iter)
        A_eff = 3 * sizeof(Float32) * (nx[iRes]+1) * (ny[iRes]) + 3 * sizeof(Float32) * (nx[iRes]) * (ny[iRes]+1) + 4 * sizeof(Float32) * (nx[iRes]) * (ny[iRes]) # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
        # t_it  = (t_toc-t_tic)/(niter-warmup_iter)   # Time per iteration [s]
        # t_it = t_toc
        T_eff = Float64(A_eff/t_it)                          # Effective memory throughput [GB/s]
        Res = Float64(nx[iRes]*ny[iRes])
        push!(T_eff_all ,  T_eff  )
        push!(Res_all   ,  Res    )
        push!(T_peak_all,  T_peak )
        push!(T_apple_all, T_apple)

        # @printf("Elapsed time = %1.3f s; No. iterations = %d; Time / iteration = %1.3f s; T_eff = %1.3f GB/s\n", t_toc-t_tic, niter-warmup_iter, t_it, T_eff/1e9)
        @printf("Nx = ny = %d; Time / iteration = %1.3e s; T_peak = %1.3f GB/s; T_eff = %1.3f GB/s\n", nx[iRes], t_it, T_peak/1e9, T_eff/1e9)
    end
    
    # Visualise and save to png
    f = Figure(resolution=(600,500))
    ax = Axis(f[1, 1], xlabel = L"nx \times ny []", ylabel = L"T [B/s]", title = "Weak scaling benchmark", limits=(minimum(Res_all), maximum(Res_all), minimum(T_eff_all), 250e9), xscale=log10, yscale=log10, xlabelsize=20, ylabelsize=20, titlesize=25, xticksize=18, yticksize=18, xminorgridvisible=true, xminorticksvisible=true, yminorgridvisible=true, yminorticksvisible=true, xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5))
    l1 = lines!( ax, Res_all, T_eff_all, linewidth=5)
    l2 = lines!( ax, Res_all, T_peak_all, linewidth=5)
    l3 = lines!( ax, Res_all, T_apple_all, linewidth=5)
    Legend(f[1,2], [l1, l2, l3], [L"T_{eff}", L"T_{peak}", L"T_{Apple}"])
    display(f)
    save("./doc/metalGPU_weakscaling.png", f, px_per_unit = 2)
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
