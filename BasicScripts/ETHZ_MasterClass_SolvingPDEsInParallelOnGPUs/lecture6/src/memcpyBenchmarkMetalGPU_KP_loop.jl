# Benchmarking array programming on Apple Silicon GPU
import Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using Metal, CairoMakie, Test, BenchmarkTools

function copy!(A, B)
    ix, iy   = thread_position_in_grid_2d()
    A[ix,iy] = B[ix,iy]
    return nothing
end

device = MTLDevice(1)
array_sizes = []
throughputs = []
for pow = 0:11
    nx = ny = 32*2^pow
    maxBuff = 3*nx*ny*sizeof(Float64)
    if (maxBuff > device.maxBufferLength) break; end
    println("Current max buffer = $(maxBuff/1e9); deviceMaxBuff = $(device.maxBufferLength/1e9)")
    A = MtlArray(zeros(Float32, nx, ny));
    B = MtlArray(rand(Float32, nx, ny));
    tic = 0.0
    iter = 1; itmax = 10000; warmup=11
    while iter<itmax
        if iter==warmup
            tic = Base.time() 
        end
        @inbounds @metal copy!(A, B)
        synchronize()
        iter += 1        
    end
    toc = Base.time()
    t_it = (toc-tic)/(itmax-warmup)
    T_tot = 2*1/1e9*nx*ny*sizeof(Float32)/t_it
    push!(array_sizes, nx)
    push!(throughputs, T_tot)
    # @printf("T_tot = %.4f\n", T_Tot)
    println("(nx=ny=$nx) T_tot = $(T_tot)")
end