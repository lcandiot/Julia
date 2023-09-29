using Plots, Plots.Measures, Printf, Pkg, BenchmarkTools
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
default(size=(600,500),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

function memcopy(; bench=:loop, compute=:ap)
    # Numerics 
    nx, ny  = 4096, 4096       # Array resolution 
    nt      = 2e4            # No. iterations
    warm_up = 11            # Warmup iterations
    # Initialisation 
    C       = rand(Float64, nx, ny) # Array 1
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
    A_eff = (3*8)*(nx)*(ny)         # Example for qDx: 1 read and 1 write for qDx, 1 read for Pf
    T_eff = A_eff/t_it              # Effective memory throughput [GB/s]  
    # Print to screen
    @printf("Elapsed time = %1.6f s; T_eff = %1.3f GB/s\n", t_toc, T_eff/1e9)
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
memcopy(; bench=:btool, compute=:ap)