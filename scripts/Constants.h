// Constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

const double dt = 0.001;
const double Chi = 0.193183325037836;
const double InitialVelocity = 0.000001;

// for elastic collisitions
//const double SpringConstant = 50;
//const double kB = 1.380649e-23; // Boltzmann constant in J/K
const double kB = 1; // Boltzmann constant in J/K

// Lennard-Jones
const double epsilon = 1;
const double sigma = 1;
const double cutoff = 2.5;

const double T_desired = 0.1;

const double Lx = 50, Ly = 50;
const int Nx = 8, Ny = 8, N = Nx * Ny;

// For Omelyan PEFRL
const double Zeta   = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi     = -0.06626458266981849;

#endif