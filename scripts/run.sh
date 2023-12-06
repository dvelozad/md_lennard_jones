#!/bin/bash

temperature=$1

# Number of particless
N_=$2
RHO_=$3

# Nx_=$2
# Ny_=$3
# Nz_=$4

# # Compile the source files
# g++ -DT=$temperature -DNX=$Nx_ -DNY=$Ny_ -DNZ=$Nz_ -c Particle.cpp
# g++ -DT=$temperature -DNX=$Nx_ -DNY=$Ny_ -DNZ=$Nz_ -fopenmp -c Collider.cpp
# g++ -DT=$temperature -DNX=$Nx_ -DNY=$Ny_ -DNZ=$Nz_ -fopenmp -c main.cpp

# # Link the object files and create the executable
# g++ -DT=$temperature -DNX=$Nx_ -DNY=$Ny_ -DNZ=$Nz_ -fopenmp Particle.o Collider.o main.o -o simul_exec_$temperature

# Compile the source files
g++ -DT=$temperature -DN_=$N_ -DRHO_=$RHO_ -c Particle.cpp
g++ -DT=$temperature -DN_=$N_ -DRHO_=$RHO_ -fopenmp -c Collider.cpp
g++ -DT=$temperature -DN_=$N_ -DRHO_=$RHO_ -fopenmp -c main.cpp

# Link the object files and create the executable
g++ -DT=$temperature -DN_=$N_ -DRHO_=$RHO_ -fopenmp Particle.o Collider.o main.o -o simul_exec_$temperature

./simul_exec_$temperature
