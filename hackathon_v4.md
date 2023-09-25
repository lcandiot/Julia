# Hackathon v4.0 project goals

During the Hackathon, my primary objective was to start developping geodynamic modeling applications using the Julia programming language. 
To achieve this, I dedicated my efforts to the development of partial differential equation (PDE) solvers tailored to purely thermal and hydrothermal problems, employing an assortment of discretization techniques.
# Model development
Throughout the week, I created models to forecast how thermal and hydrothermal systems behave in both one-dimensional (1D) and two-dimensional (2D) scenarios. 
These models encompassed solutions for the steady-state and changing conditions described by the heat equation and Darcy's law. I incorporated various time integration methods, both implicit and explicit, alongside diverse boundary conditions. 
The video below demonstrates the solution of 2D porous convection using an implicit time integration approach for temperature and fluid pressure, in addition to flux boundary conditions.

https://github.com/lcandiot/Julia/assets/50524459/f5d31b19-4bfc-4733-b62e-a20e654e82c0

# Visualisation routines
In addition to solver development, I further programmed a routine to generate automatic visual output from HPC simulations on cluster.
This program is also written in Julia and can be executed after the simulation has completed to generate mp4 files from data stored as hdf5 files.
# Outcome
The Hackathon provided an excellent opportunity for me to make rapid progress in my Julia-based geodynamic modeling projects. 
These advancements will serve as the foundation for my ongoing research endeavors. 
In essence, the event accelerated my work and set the stage for my current research initiatives. 


