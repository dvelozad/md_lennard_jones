// Collider.cpp
#include "Collider.h"
#include "Constants.h"
#include <cmath>

extern const double SpringConstant;
extern const int N;

// Collider methods
void Collider::Init(void) {
    potentialEnergy = 0;
}

void Collider::CalculateForces(Particle *particles) {
    for (int i = 0; i < N; i++) { //N + 4
        particles[i].forceX = 0; particles[i].forceY = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) { //N + 4
            Collide(particles[i], particles[j]);
        }
    }
}

void Collider::Collide(Particle &particle1, Particle &particle2) {
    double dx = particle1.x - particle2.x;
    double dy = particle1.y - particle2.y;
    double distance = sqrt(dx * dx + dy * dy);
    double overlap = particle1.radius + particle2.radius - distance;

    /*
    if (overlap > 0) {
        double forceNormal = SpringConstant * pow(overlap, 1.5);
        double forceX = forceNormal * dx / distance;
        double forceY = forceNormal * dy / distance;

        particle1.forceX += forceX;
        particle1.forceY += forceY;
        particle2.forceX -= forceX;
        particle2.forceY -= forceY;

        potentialEnergy += SpringConstant / 2.5 * pow(overlap, 2.5);
    }
    */

    double forceNormal = -24*epsilon*((pow(sigma, 6) / pow(distance, 7)) - 2*(pow(sigma, 12) / pow(distance, 13))); 
    double forceX = forceNormal * dx / distance;
    double forceY = forceNormal * dy / distance;

    particle1.forceX += forceX;
    particle1.forceY += forceY;
    particle2.forceX -= forceX;
    particle2.forceY -= forceY;
    potentialEnergy += 4*epsilon*(pow((sigma/distance), 12) - pow((sigma/distance), 6));
}

double Collider::GetPotentialEnergy(void){return potentialEnergy;};