# 2D Acoustic Wave Simulator (Finite Difference Method)

A high-performance MATLAB implementation for simulating acoustic wave propagation in 2D media using the **Finite Difference Method (FDM)**.

## Overview
This simulator solves the 2D acoustic wave equation, allowing users to define complex velocity models, various source types, and receiver arrays. It provides real-time visualization (snapshots) of the wavefield and records seismic traces (Shot Gathers) for further analysis.

## Key Features
* Numerical Core: Solves the wave equation with an optimized FDM engine (`fde.m`).
* Stability Control: Automatically calculates the time-step ($dt$) based on the CFL condition to ensure numerical stability.
* Absorbing Boundaries: Implements Absorbing Boundary Conditions (ABC) to eliminate artificial reflections from the edges.
* Source Options: Supports Ricker wavelets, sinusoidal signals, and random noise sources.
* Data Acquisition: Virtual receivers can be placed anywhere in the model to record pressure fields over time.

##  Project Structure
* `acu2Dpro.m`: The main simulation function (The Engine).
* `fde.m`: The finite difference kernel and boundary condition handler.
* `B01_exercise.m`: Layered velocity model (Reflection/Transmission).
* `B02_exercise.m`: Energy decay and geometric spreading analysis.
* `B03_exercise.m`: Beamforming and acoustic focusing on targets.
* `B04_exercise.m`: Refraction modeling and Shot Gather generation.
* `B05_exercise.m`: Complex scattering from high-velocity obstacles.
