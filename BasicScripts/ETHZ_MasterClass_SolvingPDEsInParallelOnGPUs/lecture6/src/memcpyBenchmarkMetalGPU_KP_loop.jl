# Benchmarking kernel programming on Apple Silicon GPU using a manual approach for time keeping (w/ loops).
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

# Define iteration function
function run_iterations!(A, B, iters, threads, groups)
    tic = time_ns()
    for _ = 1:iters
        Metal.@sync @metal threads=threads groups=groups memcopy!(A,B)
    end
    t_it = (time_ns() - tic)/iters * 1e-9   # * 1e-9 because t_it would be in nanoseconds
    return t_it
end

# Main function
function run_benchmark_KP(dtype=Float32)
    # Setup device and arrays
    device      = MTLDevice(1)
    array_sizes = []
    throughputs = []

    # Benchmark loop
    warm_up, iterations = 5, 35
    for pow = 0:11
        # Current resolution
        nx = ny = 32*2^pow
        maxBuff = 3*nx*ny*sizeof(dtype)

        # Sanity check, memory allocation, and initialisation
        if (maxBuff > device.maxBufferLength) break; end
        A = MtlArray(zeros(dtype, nx, ny));
        B = MtlArray(rand(dtype, nx, ny));

        # Define threads per group and no. groups
        threads = (16,16)
        groups  = cld.(size(A), threads)

        # Do warm-up iterations
        run_iterations!(A, B, warm_up, threads, groups)
        
        # Take time
        t_it = run_iterations!(A, B, (iterations-warm_up), threads, groups)
        
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