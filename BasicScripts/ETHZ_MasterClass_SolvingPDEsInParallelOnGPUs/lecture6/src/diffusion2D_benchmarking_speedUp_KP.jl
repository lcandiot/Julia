# Benchmarking 2D diffusion equation on Apple GPU
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
using IJulia
using Metal
using BenchmarkTools
using CairoMakie

# Macro definitions
@inbounds @views macro d_xa(A) esc(:( ($A[2:end  , :     ] .- $A[1:end-1, :     ]) )) end
@inbounds @views macro d_xi(A) esc(:( ($A[2:end  ,2:end-1] .- $A[1:end-1,2:end-1]) )) end
@inbounds @views macro d_ya(A) esc(:( ($A[ :     ,2:end  ] .- $A[ :     ,1:end-1]) )) end
@inbounds @views macro d_yi(A) esc(:( ($A[2:end-1,2:end  ] .- $A[2:end-1,1:end-1]) )) end
@inbounds @views macro  inn(A) esc(:( $A[2:end-1,2:end-1]                          )) end

# Define diffusion function
function diffusion2D()
    # Physics
    # lam      = Float32(1.0)                                          # Thermal conductivity
    c0       = Float32(2.0)                                          # Heat capacity
    lx, ly   = Float32(10.0), Float32(10.0)                                   # Length of computational domain in dimension x and y

    # Numerics
    nx, ny   = Int32(32768), Int32(32768)                                   # Number of gridpoints in dimensions x and y
    nt       = Int32(200)                                          # Number of time steps
    dx       = Float32(lx/(nx-1))                                    # Space step in x-dimension
    dy       = Float32(ly/(ny-1))                                    # Space step in y-dimension
    # _dx, _dy = Float32(1.0/dx), Float32(1.0/dy)

    # Setup device and arrays
    device   = MTLDevice(1)

    # Array initializations
    T    = MtlArray(zeros(Float32, nx, ny))                      # Temperature
    T2   = MtlArray(zeros(Float32, nx, ny))                      # Temperature alias
    Ci   = MtlArray(zeros(Float32, nx, ny))                      # 1/Heat capacity
    lam = _dx = _dy = dt = Float32(rand());

    # Initial conditions
    Ci .= Float32(1.0/(c0))                                              # 1/Heat capacity (could vary in space)
    T  .= MtlArray([Float32(10.0)*exp(-(((ix-1)*dx-lx/2)/2)^2-(((iy-1)*dy-ly/2)/2)^2) for ix=1:size(T,1), iy=1:size(T,2)]) # Gaussian temperature anomaly
    T2 .= T;                                                 # Assign also T2 to get correct boundary conditions.

    # Time loop
    # dt  = Float32(min(dx^2,dy^2)/lam/maximum(Ci)/4.1)                # Time step for 2D Heat diffusion
    # opts = (aspect_ratio=1, xlims=(1, nx), ylims=(1, ny), clims=(0.0, 10.0), c=:davos, xlabel="Lx", ylabel="Ly") # plotting options
    # for it = 1:nt
    #     diffusion2D_step!(T2, T, Ci, lam, dt, _dx, _dy)     # Diffusion time step.
    #     if it % 10 == 0
    #         # IJulia.clear_output(true)
    #         display(heatmap(Array(T)'; opts...))            # Visualization
    #         sleep(0.1)
    #     end
    #     T, T2 = T2, T                                       # Swap the aliases T and T2 (does not perform any array copy)
    # end

    # Benchmark
    t_it = @belapsed begin diffusion2D_step!($T2, $T, $Ci, $lam, $dt, $_dx, $_dy); end
    T_tot_lb = 3.0 /(2^30) * nx * ny * sizeof(Float32) / t_it
    println("nx*ny = $(nx*ny); T_tot_lb = $(T_tot_lb) GB/s; t_it = $(t_it) s")

    # Return
    return
end

# Define diffusion step function
function diffusion2D_step!(T2, T, Ci, lam, dt, _dx, _dy)
    threads = (16,16) # 16 by 16 seems to be the most performant configuration
    groups  = cld.(size(T2), threads)
    Metal.@sync @metal threads=threads groups=groups updateTemperature!(T2, T, Ci, lam, dt, _dx, _dy)
end
function updateTemperature!(T2, T, Ci, lam, dt, _dx, _dy)
    ix, iy = thread_position_in_grid_2d()
    nx, ny = size(T2)
    if (ix>=2 && ix<=nx-1 && iy>=2 && iy<=ny-1)
        @inbounds T2[ix,iy] = T[ix,iy] + dt * Ci[ix,iy] *
        ( 
            -( lam*_dx*( (-(T[ix+1,iy  ]-T[ix,iy])) - (-(T[ix,iy]-T[ix-1,iy  ])) ) )*_dx 
            -( lam*_dy*( (-(T[ix  ,iy+1]-T[ix,iy])) - (-(T[ix,iy]-T[ix  ,iy-1])) ) )*_dy
        )                               # Update of temperature             T_new = T_old + ∂t ∂T/∂t
    end
    return nothing
end

# Run diffusion
diffusion2D()