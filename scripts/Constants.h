// Constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Mod Verlet definitons
//const double Chi = 0.193183325037836;

// for elastic collisitions
//const double SpringConstant = 50;

// Lennard-Jones
const double kB         = 1; //1.380649e-23;
const double epsilon    = 1;
const double sigma      = 1;
const double cutoff     = 2.5;

// For Omelyan PEFRL
const double Zeta   = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi     = -0.06626458266981849;

const double Gamma      = 1;
const double T_desired  = 2;
const double Q          = 1; 

// Particles definitions
const double defaultMass = 1;
const double InitialVelocity = 0.000001;
const double Lx = 10, Ly = 10, Lz = 10;
const int Nx = 10, Ny = 10, Nz = 10, N = Nx * Ny * Nz;

// Write definitions
const double dt = 0.001;
const int timeFrame = 2;
const double totalTime = 1000;

// Steepest Descent Energy Minimization
const double minimizationStepSize = 0.0001;  
const int minimizationSteps = 0;      

// Radial distribution function definitions
const double maxDistance = 10; 
const int numBins = 200; 

#endif