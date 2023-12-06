// Constants.h

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef T
#define T 1
#endif

/*#ifndef NX
#define NX 1
#endif

#ifndef NY
#define NY 1
#endif

#ifndef NZ
#define NZ 1
#endif*/

#ifndef N_
#define N_ 1
#endif

#ifndef RHO_
#define RHO_ 1
#endif

// Simulation label
#include <cmath>
#include <string>
const std::string simulationLabel = "T" + std::string(TOSTRING(T)) + "_N" + std::string(TOSTRING(N_))+ "_RHO" + std::string(TOSTRING(RHO_));

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

const double Gamma      = 0.5;
const double T_desired  = T;
const double Q          = 1; 

// Particles definitions
const double defaultMass = 1;
const double InitialVelocity = pow(12, 0.5);
//const int Nx = NX, Ny = NY, Nz = NZ, N = Nx * Ny * Nz;
const int N = 4*N_*N_*N_;
const double L = pow(defaultMass*N / RHO_, 1. / 3.);
const double Lx = L, Ly = L, Lz = L;

// Write definitions
const double dt = 0.001;
const int timeFrame = 10;
const double totalTime = 1000;

// Steepest Descent Energy Minimization
const double minimizationStepSize = 0.0000001;  
const int minimizationSteps = 0;      

// Radial distribution function definitions
const double maxDistance = 10; 
const int numBins = 200; 

#endif