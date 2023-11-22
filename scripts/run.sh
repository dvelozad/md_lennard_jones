#!/bin/bash

# Compile the source files
g++ -c Particle.cpp
g++ -fopenmp -c Collider.cpp
g++ -fopenmp -c main.cpp

# Link the object files and create the executable
g++ -fopenmp Particle.o Collider.o main.o -o simul_exec

# Run the simulation and pipe its output to gnuplot
./simul_exec #| gnuplot