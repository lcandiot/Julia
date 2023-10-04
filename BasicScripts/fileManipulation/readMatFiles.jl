using MAT, CairoMakie, Parameters
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
# Read the MAT-file
vars = matread("/Users/lcandiot/Developer/Julia/BasicScripts/data/dummyData.mat")
@unpack d, e = vars
a = vec(d)
b = vec(e)
x = range(0, 10, length=100)
# Plotting with Makie
fig1 = Figure()
ax   = Axis(fig1[1, 1]) 
lines!(ax, a, b)
fig1
