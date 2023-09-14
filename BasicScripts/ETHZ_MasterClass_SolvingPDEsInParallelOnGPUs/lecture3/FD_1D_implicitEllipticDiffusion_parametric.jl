# Solving 1D diffusion equation to steady state with central finite differences
using GLMakie
GLMakie.activate!()
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
# Define Function
@views function steadyDiffusion_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    dc      = 1.0                   # Diffusion coefficient
    # Numerics
    ncx     = 200                   # Number of cells in x
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 100ncx                # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)     # Convergence check frequency
    fact    = 0.5:0.1:1.5           # Testing factors
    conv    = zeros(size(fact))     # Storage array for no. iterations
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    nt      = 5ncx #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    for ifact in eachindex(fact)
        # Initialisation
        re      = 2π*fact[ifact]    # Reynolds number
        ρ       = (lx/(dc*re))^2    # Numerical density
        dτ      = dx/sqrt(1/ρ) #dx^2/dc/2.1           # Time step size
        C       = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; C_ini = copy(C)
        qx      = zeros(Float64, ncx-1)
        iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
        fig1    = Figure()                 # Plotting
        ax1     = Axis(fig1[1, 1], xlabel=L"\textit{x}", ylabel=L"\textit{C}")
        ax2     = Axis(fig1[2, 1], xlabel=L"iter/ncx", ylabel=L"\textit{err}", yscale = log10)
        lines!(ax1, xc, C, color = :blue, label=L"\textit{C}")
        lines!(ax1, xc, C_ini, color = :orange, label=L"\textit{C}_{ini}")
        lines!(ax2, iter_evo, err_evo, color = :blue)
        axislegend(ax1, merge=true, unique=true)
        display(fig1) 
        # Iteration loop
        while err >= ϵtol && iter <= maxiter
            # Computation
            #qx         .-=   d\tau./(ρ.*dc + d\tau).*(qx + dc.*diff(C)./dx)                  # Ludo's formulation
            qx          .=   1.0./(dτ + ρ.*dc) .* (ρ.*dc.*qx - dc.*dτ.*diff(C)./dx) # My formulation 
            C[2:end-1] .-=   dτ.* diff(qx)./dx
            # Visualisation
            if iter % ncheck == 0
                err = maximum(abs.(diff(dc.*diff(C)./dx)./dx))
                push!(iter_evo, iter/ncx); push!(err_evo, err)
                sleep(0.1)
                empty!(ax1)
                empty!(ax2)
                lines!(ax1, xc, C, color = :blue, label=L"\textit{C}")
                lines!(ax1, xc, C_ini, color = :orange, label=L"\textit{C}_{ini}")
                lines!(ax2, iter_evo, err_evo, color = :blue)
                display(fig1)
            end
            iter += 1
        end
        conv[ifact] = iter/ncx
        # Save figure
        if ifact == length(fact)
            save("./ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/png/steadyStateDiffusion_implicit_1D.png", fig1)
        end
    end
    # Visualise result of parametric study
    fig2 = Figure()
    ax3  = Axis(fig2[1,1], xlabel="factor", ylabel="iter/ncx")
    lines!(ax3, fact, conv)  
    display(GLMakie.Screen(), fig2)
    save("./ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/png/parametricStudy.png", fig2)    
end

# Call function
steadyDiffusion_implicit_1D()