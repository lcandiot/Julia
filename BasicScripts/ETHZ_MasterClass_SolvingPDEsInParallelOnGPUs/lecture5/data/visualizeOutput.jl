if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
using CairoMakie, Printf, HDF5

# Data location
datPath = "./data"
outName = "/memcpy_kp.h5"
# Load data
nx    = h5read("$(datPath)$(outName)","/model/nx")
ny    = h5read("$(datPath)$(outName)","/model/ny")
T_eff = h5read("$(datPath)$(outName)","/monitor/T_eff")
println(T_eff)
println(nx)
println(ny)
# Visualize
fig2 = Figure(resolution=(600,500), fontsize=25)
ax1  = Axis(fig2[1,1], xlabel=L"nx \times ny", ylabel=L"T_{eff}")
lines!(ax1, nx.*ny ,T_eff/1e9)
lines!(ax1, nx.*ny, 123*ones(Float64,length(nx)))
display(fig2)