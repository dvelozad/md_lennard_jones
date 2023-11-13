// Particle.cpp
#include "Particle.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

using namespace std;

extern const double Chi;

// Particle methods


void Particle::Init(double x0, double y0, double velocityX0, double velocityY0, double mass0, double radius0) {
    x = x0; y = y0; velocityX = velocityX0; velocityY = velocityY0; mass = mass0; radius = radius0;
}

void Particle::Move_r1(double dt) {
    x += velocityX * Chi * dt;  y += velocityY * Chi * dt;
}

void Particle::Move_V(double dt) {
    velocityX += forceX * (dt / (2 * mass));  velocityY += forceY * (dt / (2 * mass));
}

void Particle::Move_r2(double dt) {
    x += velocityX * ((1 - 2 * Chi) * dt);  y += velocityY * ((1 - 2 * Chi) * dt);
}

void Particle::Draw(void) {
    cout << " , " << x << "+" << radius << "*cos(t)," << y << "+" << radius << "*sin(t)";
}

double Particle::GetX(void) { return x; }
double Particle::GetY(void) { return y; }
double Particle::GetKineticEnergy(void) { return mass * (velocityX * velocityX + velocityY * velocityY); }
double Particle::GetVelocity(void) { return sqrt(velocityX * velocityX + velocityY * velocityY); }
double Particle::GetVelocityX(void) { return velocityX; }