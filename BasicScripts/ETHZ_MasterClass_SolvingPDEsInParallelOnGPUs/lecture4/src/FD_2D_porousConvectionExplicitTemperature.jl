# Solving 1D diffusion equation to steady state with central finite differences
using CairoMakie
# Define Function
@views function steadyDiffusion_implicit_1D()
    # Physics
    lx      = 20.0                  # Length in x
    k_ηf    = 1.0                   # Diffusion coefficient
    # Numerics
    re      = 2π
    CFL     = 1.0/sqrt(2.1)
    ncx     = 200                   # Number of cells in x
    ϵtol    = 1e-8                  # Residual tolerance
    maxiter = 20ncx                  # Maximum no. iterations
    ncheck  = ceil(Int,0.25ncx)      # Convergence check frequency
    # Derived Numerics
    dx      = lx/ncx                # Spatial step size
    θ_dτ    = lx/re/CFL/dx
    β_dτ    = (re*k_ηf)/(CFL*dx*lx)
    nt      = 5ncx #ncx^2/5              # No. of time steps to compute
    xc      = LinRange(dx/2, lx-dx/2, ncx) # Coordinate array
    # Initialisation
    Pf      = @. 1.0+exp(-(xc-lx/4)^2) - xc/lx; Pf_ini = copy(Pf)
    qDx     = zeros(Float64, ncx+1)
    r_Pf    = zeros(Float64, ncx  )
    fig1    = Figure()                 # Plotting
    ax1     = Axis(fig1[1, 1])
    ax2     = Axis(fig1[2, 1], yscale = log10)
    # Iteration loop
    iter = 1; err = 2ϵtol; iter_evo = Float64[]; err_evo = Float64[]
    while err >= ϵtol && iter <= maxiter
        # Computation
        qDx[2:end-1].-=   1.0./(1.0 + θ_dτ).*(qDx[2:end-1] + k_ηf.*diff(Pf)./dx)                  # Ludo's formulation 
        Pf          .-=   diff(qDx)./dx./β_dτ
        # Visualisation
        if iter % ncheck == 0
            r_Pf .= diff(qDx)./dx
            err = maximum(abs.(r_Pf))
            push!(iter_evo, iter/ncx); push!(err_evo, err)
            sleep(0.1)
            empty!(ax1)
            empty!(ax2)
            lines!(ax1, xc, Pf, color = :blue)
            lines!(ax1, xc, Pf_ini, color = :orange)
            lines!(ax2, iter_evo, err_evo, color = :blue)
            display(fig1)
        end
        iter += 1
    end        
end

# Call function
steadyDiffusion_implicit_1D()