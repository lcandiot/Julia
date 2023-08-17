# Car travel 1D
using Plots

@views function car_travel_1D()
    # Physics
    velocity    = 113                       # Speed of the car [km/h]
    turnAround  = 200                       # Point of U-turn [km]
    travelTime  = 16                        # Travel time [h]
    direction   = 1                         # Direction of movement; + or -
    # Numerics
    deltaTime   = 0.1                               # Time step [h]
    noTimeSteps = Int(cld(travelTime,deltaTime))    # Total number of time steps

    # Array initialisation
    Time        = Array{Float64}(undef, noTimeSteps)
    Position    = Array{Float64}(undef, noTimeSteps)

    # Time loop
    for iTime = 2:noTimeSteps
        # Update time and Position 
        Time[iTime] = Time[iTime-1] + deltaTime
        Position[iTime] = Position[iTime-1] + direction*velocity*deltaTime

        # Turn arounds
        if Position[iTime] > turnAround || Position[iTime] < 0.0
            direction = -direction
        end
    end

    # Visualisation
    display(scatter(Time, Position, markersize=5,
                    xlabel="time, hrs", ylabel="distance, km",
                    framestyle=:box, legend=:none))
    return
end

car_travel_1D()