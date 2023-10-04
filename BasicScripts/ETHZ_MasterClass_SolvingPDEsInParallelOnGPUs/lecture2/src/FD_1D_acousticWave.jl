# Solving 1D transient acoustic wave equation with central finite differences
using Pkg, CairoMakie
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
# Theme
myTheme = Theme(fontsize = 25)
set_theme!(myTheme)
@views function acousticWave1D()
    # Physics
    lx      = 20.0
    ρ, β    = 1.0, 1.0

    # Numerics
    ncx     = 200                   # Number of cells in x
    nVis    = 2                     # Visualise every XXX
    
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    dt      = dx/sqrt(1/β/ρ)        # Time step size
    nt      = ncx^2/100             # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    Pa      = @. exp(-(xc-lx/4)^2); Pa_ini = copy(Pa) # Acoustic pressure
    vx      = zeros(Float64, ncx-1)
    fig1    = Figure()                 # Plotting
    ax      = Axis(fig1[1, 1], xlabel=L"\textit{x} [m]", ylabel=L"\textit{P}_{a} [Pa]", limits=(0, lx, -1, 1), title="1D Acoustic Wave")
    lines!(ax, xc, Pa)
    lines!(ax, xc, Pa_ini)
    display(fig1) 
    # Time loop
    for it = 1:nt
        # Computation
        vx          .-= dt./ρ.*diff(Pa)./dx
        Pa[2:end-1] .-= dt./β.*diff(vx)./dx
        # Visualisation
        if it%nVis == 0
            sleep(0.05) 
            empty!(ax)
            lines!(ax, xc, Pa, color=:blue, xlim=(0, lx))
            lines!(ax, xc, Pa_ini, color=:orange, xlim=(0, lx))
            display(fig1)
        end
    end 
end

# Call function
acousticWave1D()