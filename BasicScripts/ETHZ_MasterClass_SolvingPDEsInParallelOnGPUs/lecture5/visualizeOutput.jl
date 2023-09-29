using CairoMakie, Printf, HDF5, Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

struct dataStructure
    nx::Vector{Int64}
    ny::Vector{Int64}
    T_eff::Vector{Float64}
end

function visualizeOutput(; printFig=true)
    nx = 16 * 2 .^(1:8)
    # Data location
    datPath = "./data/"
    docPath = "./doc/"
    outName = ["Pf_diffusion_2D_perf.h5", "Pf_diffusion_2D_perf_loop.h5", "Pf_diffusion_2D_perf_loop_fun.h5", "memcpy_ap.h5", "memcpy_kp.h5"]
    
    # Prepare figure
    fig2 = Figure(resolution=(600,500), fontsize=20, figure_padding=15)
    ax1  = Axis(fig2[1,1], xlabel=L"nx \times ny", ylabel=L"T_{eff}", xscale=log10)
    # Load and plot my data 
    for iFile in eachindex(outName) 
        fname = "$(datPath)$(outName[iFile])"
        dat = readThroughput(fname)
        # Visualize
        lines!(ax1, dat.nx.*dat.ny ,dat.T_eff/1e9, label=outName[iFile])
        if iFile==1
            # Add maximum loop bandwidth and vendor announced values
            lines!(ax1, dat.nx.*dat.ny , 123*ones(Float64,length(dat.nx)), linestyle=:dash, label="Max loop")
            lines!(ax1, dat.nx.*dat.ny , 200*ones(Float64,length(dat.nx)), linestyle=:dash, label="Vendor")
        end
    end
    # Add legend and display
    fig2[2,1] = Legend(fig2, ax1, "Solvers", framevisible=false, orientation=:horizontal, labelsize=12, nbanks=3)
    display(fig2)
    # Print figure
    if printFig
        save("$(docPath)scalingTest.png", px_per_unit=2, fig2)
    end
    # Return
    return
end

# Read input function
function readThroughput(fname)
    # Store data in struct
    dat = dataStructure(h5read(fname,"/model/nx"), 
    h5read(fname,"/model/ny"), 
    h5read(fname,"/monitor/T_eff"))
    # Return data struct
    return dat
end

# Call main 
visualizeOutput(; printFig=true)