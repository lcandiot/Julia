# Parallel programming benchmark 
## Description
Different implementations of a Darcy equation solver have been developed. These include an optimized Julia [array programming](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture5/src/l5_Pf_diffusion_2D_perf.jl) version in which all division operations of flux calculations and pressure update have been replaced by inverse multiplications. Based on this implementation, another version has been developed. There the flux calculations and pressure update are executed using [loops](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture5/src/l5_Pf_diffusion_2D_perf_loop.jl). Julia is already optimized for loops and thus an increase in performance compared to array programming is expected. Finally, the flux calculations and the pressure update have been wrapped into [kernels](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture5/src/l5_Pf_diffusion_2D_perf_loop_fun.jl). Multi-threading has been activated for the loop versions and simulations were performed on 10 Julia Threads. The [Julia BenchmarkTools](https://juliaci.github.io/BenchmarkTools.jl/stable/) have been used for evaluation and the benchmark has been performed on an Apple Silicon M1 Pro chip (2021). For comparison, results are evaluated against memory copying operations (both for array and kernel programming). Memory copying is measured by taking the elapsed time for computation manually. The figure below shows the benchmark results.  
![scalingTest](https://github.com/lcandiot/Julia/assets/50524459/0fbd23b4-32ed-4475-bdf2-e667006c6f11)
## Results
Of all tested solver implementations, kernel programming achieves results in the least amount of computing time and at the highest bandwidth (ca. 100 GB/s). Peak bandwidth for memory copying on this machine is predicted at 123 GB/s, which is almost only half of the vendor announced limit of 200 GB/s. Reasons for this discrepancy could be that the vendor limit is only reached for GPU applications. Further testing is required to verify this hypothesis.

# Unit testing
A unit test has been conducted using the [kernel solver](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture5/src/l5_Pf_diffusion_2D_perf_loop_fun.jl) implementation. At three different locations within the domain the pressure values are extracted and compared against reference values at different resolutions. The test is passed for all resolutions (see code output below).
```Test Summary: | Pass  Total  Time
Pf test set   |    1      1  0.0s
Test Summary: | Pass  Total  Time
Pf test set   |    1      1  0.0s
Test Summary: | Pass  Total  Time
Pf test set   |    1      1  0.0s
Test Summary: | Pass  Total  Time
Pf test set   |    1      1  0.0s
```


