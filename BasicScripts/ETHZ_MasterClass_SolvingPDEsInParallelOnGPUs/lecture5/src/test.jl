using Printf
function indexingCheck()

    nx, ny = 512, 512
    A      = rand(Float64, nx, ny)
    B      = copy(A)
    C      = copy(A)
    @time x_innerIndex(A, B, C)
    @time y_innerIndex(A, B, C)
    # Return
    return
end

function x_innerIndex(A, B, C)
    nx, ny = size(A)
    for iy=1:ny 
        for ix=1:nx 
            C[ix,iy] = A[ix,iy] + B[ix,iy]
        end
    end

    # Return
    return nothing
end

function y_innerIndex(A, B, C)
    nx, ny = size(A)
    for ix=1:nx 
        for iy=1:ny 
            C[ix,iy] = A[ix,iy] + B[ix,iy]
        end
    end

    # Return
    return nothing
end

# Call function
indexingCheck()