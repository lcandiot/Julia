# Benchmarking kernel programming on Apple Silicon GPU using Julias benchmarking tools
import Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using Metal, CairoMakie, Test, BenchmarkTools

# Define copy operation function
function memcopy!(A, B)
    ix, iy   = thread_position_in_grid_2d()
    @inbounds A[ix,iy] = B[ix,iy]
    return nothing
end

# Define run memcopy function
function run_memcopy!(A, B, threads, groups)
    Metal.@sync @metal threads=threads groups=groups memcopy!(A,B)
end

# Main function
function run_benchmark_KP(dtype=Float32)
    # Setup device and arrays
    device      = MTLDevice(1)
    array_sizes = []
    throughputs = []
    # Benchmark loop
    for pow = 0:11
        # Current resolution
        nx = ny = 32*2^pow
        # Sanity check, memory allocation, and initialisation
        maxBuff = 3*nx*ny*sizeof(dtype)
        if (maxBuff > device.maxBufferLength) break; end
        A = MtlArray(zeros(dtype, nx, ny));
        B = MtlArray(rand(dtype, nx, ny));
        # Define threads per group and no. groups
        threads = (32,32)
        groups  = cld.(size(A), threads)
        # Take time
        t_it = @belapsed begin run_memcopy!($A, $B, $threads, $groups) end
        # Calculate memory throughput
        T_tot = 2/2^30*nx*ny*sizeof(dtype)/t_it
        push!(array_sizes, nx)
        push!(throughputs, T_tot)
        println("(nx=ny=$nx) T_tot = $(T_tot)")
    end
    return
end

# Run script
run_benchmark_KP()