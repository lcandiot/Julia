# Julia

This repository contains a divers set of algorithms written in Julia language.
These perform different tasks from simple function exercises up to solving PDEs. 

# Steady-state diffusion solver

Solves the 1D diffusion equation to a steady-state using an implicit pseudo time iteration scheme ->
[Steady-state chemical diffusion](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/FD_1D_implicitEllipticDiffusion_parametric.jl)
![Alt text](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture3/png/steadyStateDiffusion_implicit_1D.png?raw=true)

# Porous convection 2D
Two codes have been developed treating the time integration either explicitly or implicitly. The temperature perturbation in the middle of the domain heated from bottom and cooled from top induces local differences in the specific gravity of the fluid. Hot fluid from the bottom rises and cold fluid from the top sinks down producing a convective pattern.-> [Porous convection 2D explicit time integration](BasicScripts/ETHZ_MasterClass_SolvingPDEsInParallelOnGPUs/lecture4/src/FD_2D_porousConvectionExplicitTemperature.jl)

https://github.com/lcandiot/Julia/assets/50524459/b88e72c3-e671-46f3-adda-a719291e9fcd

