# St. Venant Equations Solver

This repository contains code developed as part of a class project in 2011 for solving the St. Venant equations. The project includes both MATLAB and Python implementations of the solver.

## Project Overview

The St. Venant equations are shallow water equations used to model the flow of water in open channels. This project aimed to develop a numerical solver for these equations, and it includes MATLAB and Python implementations.

### Files and Directories

#### Matlab
- `matlab/`: Contains MATLAB code for the solver.
  - `make_vectors.m`: The make_vectors function is responsible for setting up the essential vectors and arrays required for the St. Venant equations solver. It prepares the spatial and temporal grids, initializes the state variables, and allocates memory for various arrays. 
  - `ArtificialViscosity.m`: The "ArtificialViscosity" MATLAB function is used within the St. Venant equations solver to incorporate artificial viscosity into the numerical simulation. Artificial viscosity is a technique employed in numerical methods for simulating fluid flow to improve the stability and accuracy of the calculations. This function computes and applies artificial viscosity to the input state variable U. 
  - `corrector.m`: The corrector function complements the predictor function by further refining the solution in the St. Venant equations solver. It also takes inputs such as the current state, time step information, and grid indices.
  - `predictor.m`: The predictor function is responsible for advancing the solution in time for the St. Venant equations solver. It takes several inputs, including the current state of the solution, time step information, and grid indices.
  - `single_wave.m`: The single_wave MATLAB script is part of the St. Venant equations solver and serves as the main driver for a simple example simulation. It sets up the initial conditions, defines parameters, and iterates over time to solve the equations.
  - `harbor.m`: This script simulates the behavior of water flow in a harbor and includes various components to handle wave generation, boundary conditions, and numerical computations.

#### Python
- `python/`: Contains Python code for the solver.
  - `st_v.py`: This code consists of a class called SAINT_VENANT that performs various operations related to the simulation.
  - `sv_config.yml`: Configuration file needed to run the solver.

- `docs/`: Contains project documentation.
  - `CEE279_Project2_FrameJonathan_20110608b.pdf`: A simple class writeup, no judgements allowed.
