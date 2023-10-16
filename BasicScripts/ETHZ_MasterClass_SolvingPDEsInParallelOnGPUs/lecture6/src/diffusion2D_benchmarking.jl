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
    # nt       = Int32(100)                                          # Number of time steps
    dx       = Float32(lx/(nx-1))                                    # Space step in x-dimension
    dy       = Float32(ly/(ny-1))                                    # Space step in y-dimension
    # _dx, _dy = Float32(1.0/dx), Float32(1.0/dy)

    # Setup device and arrays
    device   = MTLDevice(1)

    # Array initializations
    T    = MtlArray(zeros(Float32, nx, ny))                      # Temperature
    Ci   = MtlArray(zeros(Float32, nx, ny))                      # 1/Heat capacity
    qTx  = MtlArray(zeros(Float32, nx-1, ny-2))                  # Heat flux, x component
    qTy  = MtlArray(zeros(Float32, nx-2, ny-1))                  # Heat flux, y component
    dTdt = MtlArray(zeros(Float32, nx-2, ny-2))                  # Change of Temperature in time
    lam = _dx = _dy = dt = Float32(rand());

    # Initial conditions
    Ci .= Float32(1.0/(c0))                                              # 1/Heat capacity (could vary in space)
    T  .= MtlArray([Float32(10.0)*exp(-(((ix-1)*dx-lx/2)/2)^2-(((iy-1)*dy-ly/2)/2)^2) for ix=1:size(T,1), iy=1:size(T,2)]) # Initialization of Gaussian temperature anomaly

    # Time loop
    # dt  = Float32(min(dx^2,dy^2)/lam/maximum(Ci)/4.1)                # Time step for 2D Heat diffusion
    # opts = (aspect_ratio=1, xlims=(1, nx), ylims=(1, ny), clims=(0.0, 10.0), c=:davos, xlabel="Lx", ylabel="Ly") # plotting options
    # for it = 1:nt
    #     diffusion2D_step!(T, Ci, qTx, qTy, dTdt, lam, dt, _dx, _dy) # Diffusion time step.
    #     if it % 10 == 0
    #         #IJulia.clear_output(true)
    #         Plots.display(Plots.heatmap(Array(T)'; opts...))            # Visualization
    #         sleep(0.1)
    #     end
    # end

    # Benchmark
    t_it = @belapsed begin runBenchmark!($T, $Ci, $qTx, $qTy, $dTdt, $lam, $dt, $_dx, $_dy); end
    T_tot_lb = 11.0 /(2^30) * nx * ny * sizeof(Float32) / t_it
    T_eff    = 3.0 /(2^30) * nx * ny * sizeof(Float32) / t_it
    println("nx*ny = $(nx*ny); T_tot_lb = $(T_tot_lb) GB/s; T_eff = $(T_eff) GB/s; t_it = $(t_it) s")
    t_it_task1     = t_it
    T_tot_lb_task1 = T_tot_lb

    # Return
    return T_tot_lb_task1, t_it_task1
end

# Define diffusion step function
function diffusion2D_step!(T, Ci, qTx, qTy, dTdt, lam, dt, _dx, _dy)
    # qTx     .= .-lam.*@d_xi(T).*_dx                              # Fourier's law of heat conduction: qT_x  = -λ ∂T/∂x
    # qTy     .= .-lam.*@d_yi(T).*_dy                              # ...                               qT_y  = -λ ∂T/∂y
    # dTdt    .= @inn(Ci).*(.-@d_xa(qTx).*_dx .- @d_ya(qTy).*_dy)  # Conservation of energy:           ∂T/∂t = 1/cp (-∂qT_x/∂x - ∂qT_y/∂y)
    # @inn(T) .= @inn(T) .+ dt.*dTdt                               # Update of temperature             T_new = T_old + ∂t ∂T/∂t
    @inbounds @views qTx     .= .-lam.*(T[2:end,  2:end-1].-T[1:end-1,2:end-1]).*_dx                              # Fourier's law of heat conduction: qT_x  = -λ ∂T/∂x
    @inbounds @views qTy     .= .-lam.*(T[2:end-1,2:end  ].-T[2:end-1,1:end-1]).*_dy                              # ...                               qT_y  = -λ ∂T/∂y
    @inbounds @views dTdt    .= (Ci[2:end-1,2:end-1]).*(.-(qTx[2:end,:].-qTx[1:end-1,:]).*_dx .- (qTy[:,2:end].-qTy[:,1:end-1]).*_dy)  # Conservation of energy:           ∂T/∂t = 1/cp (-∂qT_x/∂x - ∂qT_y/∂y)
    @inbounds @views T[2:end-1,2:end-1] .= T[2:end-1,2:end-1] .+ dt.*dTdt                    
    return nothing
end

# Define run function
function runBenchmark!(T, Ci, qTx, qTy, dTdt, lam, dt, _dx, _dy) 
    Metal.@sync diffusion2D_step!(T, Ci, qTx, qTy, dTdt, lam, dt, _dx, _dy)
end

# Run diffusion
diffusion2D()