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

# Porous convection 2D
Two codes have been developed treating the time integration either explicitly or implicitly. The temperature perturbation in the middle of the domain heated from bottom and cooled from top induces local differences in the specific gravity of the fluid. Hot fluid from the bottom rises and cold fluid from the top sinks down producing a convective pattern.-> [Porous convection 2D explicit time integration](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/src/FD_2D_porousConvectionExplicitTemperature.jl)

https://github.com/lcandiot/Julia/assets/50524459/b88e72c3-e671-46f3-adda-a719291e9fcd

