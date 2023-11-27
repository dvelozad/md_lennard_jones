// Particle.cpp
#include "Particle.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

using namespace std;

extern const double Chi;

// Particle methods
void Particle::ApplyPeriodicBoundaryConditions(double Lx, double Ly) {
    if (x < 0) x += Lx;
    else if (x > Lx) x -= Lx;

    if (y < 0) y += Ly;
    else if (y > Ly) y -= Ly;

    /*    if (x < 0) velocityX = -1*velocityX;
    else if (x > Lx) velocityX = -1*velocityX;

    if (y < 0) velocityY = -1*velocityY;
    else if (y > Ly) velocityY = -1*velocityY;*/

}

void Particle::Init(double x0, double y0, double velocityX0, double velocityY0, double mass0, double radius0) {
    x = x0; y = y0; velocityX = velocityX0; velocityY = velocityY0; mass = mass0; radius = radius0;
}

void Particle::Move_r1(double dt, double constant) {
    x += velocityX * constant * dt;  y += velocityY * constant * dt;
    //x += velocityX * Chi * dt;  y += velocityY * Chi * dt;
    ApplyPeriodicBoundaryConditions(Lx, Ly);
}

void Particle::Move_V(double dt, double constant) {
    velocityX += (forceX / mass) * constant * dt;  velocityY += (forceY / mass) * constant * dt;
    //velocityX += (forceX / mass) * (1/2) * dt;  velocityY += (forceY / mass) * (1/2) * dt;
}

void Particle::Draw(void) {
    cout << " , " << x << "+" << radius << "*cos(t)," << y << "+" << radius << "*sin(t)";
}

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

double Particle::GetMass(void) { return mass; }
double Particle::GetX(void) { return x; }
double Particle::GetY(void) { return y; }
double Particle::GetKineticEnergy(void) { return 0.5*mass * (velocityX * velocityX + velocityY * velocityY); }
double Particle::GetVelocity(void) { return sqrt(velocityX * velocityX + velocityY * velocityY); }
double Particle::GetVelocityX(void) { return velocityX; }
double Particle::GetVelocityY(void) { return velocityY; }


void Particle::MinimizeEnergy(double minimizationStep) {
    // Calculate the magnitude of the force
    double forceMagnitude = sqrt(forceX * forceX + forceY * forceY);

    // Scale the step size based on the force magnitude
    double scaledStep = minimizationStep;
    if (forceMagnitude > 0) {
        scaledStep = minimizationStep / forceMagnitude;
    }

    // Update positions
    x -= forceX * scaledStep;
    y -= forceY * scaledStep;

    // Apply periodic boundary conditions
    ApplyPeriodicBoundaryConditions(Lx, Ly);
}

void Particle::RescaleVelocity(double lambda) {
    velocityX *= lambda;
    velocityY *= lambda;
}