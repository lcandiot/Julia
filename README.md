# Julia

This repository contains a divers set of algorithms written in the Julia language.
The developed programs showcase different types like tuples, sets or functions, partial differential equation solvers for multidimensional problems, as well as scaling benchmarks, and unit testing.

## Types
A basic example of a function that calculates the sum of numbers in Julia can be programmed as
```
function summit(args...)
    sum = 0
    for a in args
        sum += a 
    end

    return sum
end
println(summit(1,2))
```
This code snippet can be copied and executed in a Julia REPL. More examples can be found in -> [Types](BasicScripts/types).

## Calculate the acoustic pressure
This algorithm solves the acoustic wave propagation in a 1D domain -> [Acoustic wave](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture2/src/FD_1D_acousticWave.jl).

## Steady-state diffusion solver

Solving the coupled 1D reaction-diffusion equation to a steady-state using an implicit pseudo time iteration scheme ->
[Steady-state reaction diffusion](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/src/FD_1D_implicitSteadyDiffusionReaction.jl).
The image below shows the result of a steady-state diffusion simulation as well as the residual convergence.


![Alt text](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/doc/steadyStateDiffusion_implicit_1D.png?raw=true)

## Porous convection 2D
Two codes have been developed treating the time integration of temperature either [explicitly](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/src/FD_2D_porousConvectionExplicitTemperature_BConFluxes.jl) or [implicitly](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/src/FD_2D_porousConvection_implicit_BConFluxes.jl).

The movie below shows the simulation of convection in a porous medium in 2D. A thermal anomaly is emplaced in the center of the domain to kickstart convection. The system is heated from below and cooled from above to maintain convection over time. Boundary conditions are implemented on the fluxes, i.e. directly on the model boundary.

<video autoplay loop src=https://github.com/lcandiot/Julia/assets/50524459/a683f976-da68-4e8b-906a-fefc923e4d44 type="video/mp4">
</video>

