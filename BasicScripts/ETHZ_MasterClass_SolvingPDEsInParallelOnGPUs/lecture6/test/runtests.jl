using Plots, Test, HDF5

Pf_cpu = h5open("./test/cpu_testData.h5", "r") do outfile
    read(outfile,"centers/Pf")
end
Pf_gpu = h5open("./test/gpu_testData.h5", "r") do outfile
    read(outfile,"centers/Pf")
end
Plots.display(Plots.heatmap((Pf_cpu-Pf_gpu);aspect_ratio=1,c=:turbo))
# Pf_test = maximum(abs.(Pf_cpu))-maximum(abs.(Pf_gpu))

@testset "Pf_diffusion_2D_cpuVsgpu" begin
    @test maximum(abs.(Pf_cpu))-maximum(abs.(Pf_gpu)) ≈ 0.0
end


