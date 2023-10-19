# Parallel programming on GPU
This directory contains parallelised routines running on both Mac and NVIDIA GPUs. For the former, routines rely on the [Metal.jl](https://github.com/JuliaGPU/Metal.jl) package. 

## Running on Apple Silicon GPU
In order to evaluate the performance of algorithms, a memory copy benchmark is performed first. Peak bandwidth is reached at ca. 170 GB/s which is 15 % lower compared to the vendor limit.

In following a benchmark for 2D diffusion algorithms is performed. Both array and kernel programming approaches have been tested using a grid resolution of ca. 1B grid cells. A single iteration is performed and effektive memory throughput is calculated as 

```
T_eff = 2*D_u + D_k
```
whereby `D_u` is the number of unknown fields (i.e., that depend on their history) and `D_k` is the number of known fields (variables that do not depend on their history). 

Maximum effective throughput is ca. 163 GB/s and the computational time for a single iteration is ca. 0.07 s. Reaching this limit requires kernel programming.

Second, a 2D fluid pressure diffusion solver is benchmarked. This routine relies on kernel programming, and a weak scaling benchmark is performed, i.e. the resolution and, hence, maximum number of threads is increased and the speedup is measured via the effective memory throughput. The figure below shows the benchmark results.

![metalGPU_weakscaling](https://github.com/lcandiot/Julia/assets/50524459/93153270-b462-4166-8e29-7c1118ad77ee)

The fluid pressure routine reaches a higher bandwidth compared to simple memory benchmark. The maximum bandwidth reached before running out of memory is ca. 190 GB/s which is close to the maximum bandwidth announced by Apple (green line).
