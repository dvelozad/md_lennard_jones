// Particle.cpp
#include "Particle.h"
#include "Constants.h"
#include <cmath>
#include <random>
#include <iostream>


// Random number generation
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> dis(0, 1); // Standard normal distribution


using namespace std;

extern const double Chi;

// Particle methods
void Particle::ApplyPeriodicBoundaryConditions(double Lx, double Ly, double Lz) {
    if (x < 0) x += Lx;
    else if (x > Lx) x -= Lx;

    if (y < 0) y += Ly;
    else if (y > Ly) y -= Ly;

    if (z < 0) z += Lz;
    else if (z > Lz) z -= Lz;

    /*    if (x < 0) velocityX = -1*velocityX;
    else if (x > Lx) velocityX = -1*velocityX;

    if (y < 0) velocityY = -1*velocityY;
    else if (y > Ly) velocityY = -1*velocityY;*/

}

void Particle::Init(double x0, double y0, double z0, double velocityX0, double velocityY0, double velocityZ0, double mass0, double radius0) {
    x = x0; y = y0; z = z0; velocityX = velocityX0; velocityY = velocityY0; velocityZ = velocityZ0; mass = mass0; radius = radius0;
}

void Particle::Move_r1(double dt, double constant) {
    x += velocityX * constant * dt;  
    y += velocityY * constant * dt;
    z += velocityZ * constant * dt;
    //x += velocityX * Chi * dt;  y += velocityY * Chi * dt;
    ApplyPeriodicBoundaryConditions(Lx, Ly, Lz);
}

void Particle::Move_V(double dt, double constant) {
    velocityX += (forceX / mass - xi * velocityX) * constant * dt;  
    velocityY += (forceY / mass - xi * velocityY) * constant * dt;
    velocityZ += (forceZ / mass - xi * velocityZ) * constant * dt;
    //velocityX += (forceX / mass) * constant * dt;  velocityY += (forceY / mass) * constant * dt;
    //velocityX += (forceX / mass) * (1/2) * dt;  velocityY += (forceY / mass) * (1/2) * dt;
}

void Particle::Draw(void) {
    cout << " , " << x << "+" << radius << "*cos(t)," << y << "+" << radius << "*sin(t)";
}

/*
void Particle::Move_r1_(double dt) {
    x += velocityX * Chi * dt;  y += velocityY * Chi * dt;
    ApplyPeriodicBoundaryConditions(Lx, Ly);
}

void Particle::Move_r2_(double dt) {
    x += velocityX * ((1 - 2 * Chi) * dt);  y += velocityY * ((1 - 2 * Chi) * dt);
    ApplyPeriodicBoundaryConditions(Lx, Ly);
}

void Particle::Move_V_(double dt) {
    velocityX += (forceX / mass) * (1/2) * dt;  velocityY += (forceY / mass) * (1/2) * dt;
}
*/

double Particle::GetMass(void) { return mass; }
double Particle::GetX(void) { return x; }
double Particle::GetY(void) { return y; }
double Particle::GetZ(void) { return z; }
double Particle::GetKineticEnergy(void) { return 0.5*mass * (velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ); }
double Particle::GetVelocity(void) { return sqrt(velocityX * velocityX + velocityY * velocityY + velocityZ * velocityZ); }
double Particle::GetVelocityX(void) { return velocityX; }
double Particle::GetVelocityY(void) { return velocityY; }
double Particle::GetVelocityZ(void) { return velocityZ; }


void Particle::MinimizeEnergy(double minimizationStep) {
    // Calculate the magnitude of the force
    double forceMagnitude = sqrt(forceX * forceX + forceY * forceY + forceZ * forceZ);

    // Scale the step size based on the force magnitude
    double scaledStep = minimizationStep;
    if (forceMagnitude > 0) {
        scaledStep = minimizationStep / forceMagnitude;
    }

    // Update positions
    x -= forceX * scaledStep;
    y -= forceY * scaledStep;
    z -= forceZ * scaledStep;

    // Apply periodic boundary conditions
    ApplyPeriodicBoundaryConditions(Lx, Ly, Lz);
}

void Particle::RescaleVelocity(double lambda) {
    velocityX *= lambda;
    velocityY *= lambda;
    velocityZ *= lambda;
}

void Particle::UpdateThermostat(double kineticEnergy){
    xidot += (kineticEnergy - (3/2) * N * kB * T_desired) * dt / Q;
    xi += xidot * dt;
}

void Particle::UpdateVelocity(double dt, double gamma, double temperature) {

    // Calculate the random force magnitude
    double randomForceMagnitude = sqrt(2.0 * this->mass * kB * temperature * gamma * dt);
    
    // Update velocity
    this->velocityX += (this->forceX / this->mass - gamma * this->velocityX) * dt 
                       + randomForceMagnitude * dis(gen) / this->mass;
    this->velocityY += (this->forceY / this->mass - gamma * this->velocityY) * dt 
                       + randomForceMagnitude * dis(gen) / this->mass;
    this->velocityZ += (this->forceZ / this->mass - gamma * this->velocityZ) * dt 
                       + randomForceMagnitude * dis(gen) / this->mass;
}

double Particle::DistanceTo(Particle &other){
    double dx = x - other.x;    
    double dy = y - other.y;
    double dz = z - other.z;

    dx -= Lx * round(dx / Lx);
    dy -= Ly * round(dy / Ly);
    dz -= Lz * round(dz / Lz);

    double distance = sqrt(dx * dx + dy * dy + dz * dz);

    return distance;
}