using MAT
using GLMakie
using Parameters

# Read the MAT-file
vars = matread("./data/dummyData.mat")
@unpack d, e = vars
a = vec(d)
b = vec(e)
x = range(0, 10, length=100)
# Plotting with Makie
fig1 = Figure()
ax   = Axis(fig1[1, 1]) 
lines!(ax, a, b)
fig1
