// Constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

const double dt = 0.0001;
const double Chi = 0.193183325037836;
const double InitialVelocity = 0.001;

// for elastic collisitions
//const double SpringConstant = 50;

// Lennard-Jones
const double epsilon = 0.0103;
const double sigma = 3.4;
const double cutoff = 2.5*sigma;

const double Lx = 100, Ly = 100;
const int Nx = 6, Ny = 6, N = Nx * Ny;

// For Omelyan PEFRL
const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

#endif
