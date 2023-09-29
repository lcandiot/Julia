using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools, CairoMakie, HDF5
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
#default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
fig1 = Figure(resolution=(600,500), fontsize=25)
ax1  = Axis(fig1[1,1], xlabel=L"nx \times ny", ylabel=L"T_{eff}") 
function memcopy(; bench=:loop, compute=:kp, printFig=true, writeOut=true)
    # Numerics 
    nx = ny  = 16 * 2 .^(1:8)               # Array resolution 
    nt      = 5e2                           # No. iterations
    warm_up = 11                            # Warmup iterations
    T_eff = zeros(Float64, length(nx))
    for iRes in eachindex(nx)
        # Initialisation 
        C       = rand(Float64, nx[iRes], ny[iRes]) # Array 1
        C2      = copy(C)               # Array 2
        A       = copy(C)               # Array 3
        A_eff   = 0.0
        t_it    = 0.0
        if bench == :loop
            # Start iteration loop
            t_tic = 0.0
            for it=1:nt
                if it==warm_up
                    t_tic = Base.time()
                end
                if compute == :ap
                    compute_AP!(C, C2, A)   
                elseif compute == :kp 
                    compute_KP!(C, C2, A)             
                end
            end
            t_toc = Base.time() - t_tic
            t_it  = t_toc/(nt-warm_up)
        elseif bench == :btool
            # Use benchmark tool
            if compute == :ap
                t_toc = @belapsed compute_AP!($C, $C2, $A)
            elseif compute == :kp
                t_toc = @belapsed compute_KP!($C, $C2, $A)
            end
            t_it = t_toc
        end
        # Monitor quantities
        A_eff       = (3*8)*(nx[iRes])*(ny[iRes])         # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
        T_eff[iRes] = A_eff/t_it              # Effective memory throughput [GB/s]  
        # Print to screen
        @printf("Elapsed time = %1.6f s; T_eff = %1.3f GB/s\n", t_toc, T_eff[iRes]/1e9)
    end
    # Visualise
    lines!(ax1, nx.*ny ,T_eff/1e9)
    if compute==:kp
        lines!(ax1, nx.*ny, 123*ones(Float64,length(nx)))
    end
    display(fig1)
    # Print figure
    if printFig && compute == :ap
        save("./doc/memcpScaling_ap.png",px_per_inch=2, fig1) 
    elseif printFig && compute == :kp
        save("./doc/memcpScaling_kp.png",px_per_inch=2, fig1)        
    end
    # Write output
    if writeOut && compute == :ap
        h5write("./data/memcpy_ap.h5", "monitor/T_eff", T_eff)
        h5write("./data/memcpy_ap.h5", "model/nx",      nx   )
        h5write("./data/memcpy_ap.h5", "model/ny",      ny   ) 
    elseif writeOut && compute == :kp
        h5write("./data/memcpy_kp.h5", "monitor/T_eff", T_eff)
        h5write("./data/memcpy_kp.h5", "model/nx",      nx   )
        h5write("./data/memcpy_kp.h5", "model/ny",      ny   )          
    end
    # Return
    return 
end

# Array programming 
function compute_AP!(C, C2, A)
    C2 .= C .+ A
    # Return
    return nothing
end

# Kernel programming
function compute_KP!(C, C2, A)
    nx ,ny = size(C)
    Threads.@threads for iy=1:ny
        for ix = 1:nx
            C2[ix,iy] = C[ix,iy] + A[ix,iy]
        end
    end
    # Return
    return nothing
end

# Call main
memcopy(; bench=:btool, compute=:kp, printFig=false, writeOut=true)