# Benchmarking array programming on Apple Silicon GPU
import Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using Metal, CairoMakie, Test, BenchmarkTools

function copy!(A, B)
    ix, iy   = thread_position_in_grid_2d()
    A[ix,iy] = B[ix,iy]
    return
end

device = MTLDevice(1)
array_sizes = []
throughputs = []
for pow = 0:11
    nx = ny = 32*2^pow
    if (3*nx*ny*sizeof(Float64) > Float64(device.maxBufferLength)) break; end
    A = MtlArray(zeros(Float32, nx, ny));
    B = MtlArray(rand(Float32, nx, ny));
    t_it = @belapsed begin 
        @metal threads=nx*ny copy!($A, $B)
        synchronize() 
    end
    T_tot = 2*1/1e9*nx*ny*sizeof(Float32)/t_it
    push!(array_sizes, nx)
    push!(throughputs, T_tot)
    println("(nx=ny=$nx) T_tot = $(T_tot)")
end