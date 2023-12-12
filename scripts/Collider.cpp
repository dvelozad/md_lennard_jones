// Collider.cpp
#include "Collider.h"
#include "Constants.h"
#include <cmath>
#include <omp.h>
#include <iostream>

extern const double SpringConstant;
extern const int N;

// Collider methods
void Collider::Init(void) {
    potentialEnergy = 0;
}

void Collider::CalculateForces(Particle *particles) {
    potentialEnergy = 0;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) { //N + 4
        particles[i].forceOX = particles[i].forceX; particles[i].forceOY = particles[i].forceY; particles[i].forceOZ = particles[i].forceZ;
        particles[i].forceX = 0; particles[i].forceY = 0; particles[i].forceZ = 0;
    }

    /*    
    #pragma omp parallel for 
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) { //N + 4
            Collide(particles[i], particles[j]);
        }
    }*/

    // Update neighbor lists if necessary
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        particles[i].UpdateNeighborList(particles);
    }

    // Calculate forces
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int j : particles[i].neighborList) {
            if (i != j) { 
                Collide(particles[i], particles[j]);
            }
        }
    }

}

void Collider::Collide(Particle &particle1, Particle &particle2) {
    double dx = particle1.x - particle2.x;
    double dy = particle1.y - particle2.y;
    double dz = particle1.z - particle2.z;

    dx -= Lx * round(dx / Lx);
    dy -= Ly * round(dy / Ly);
    dz -= Lz * round(dz / Lz);

    double distance = sqrt(dx * dx + dy * dy + dz * dz);
   
   if (distance < cutoff) { 
        // Force at cutoff
        double forceCutoff = -24*epsilon*((pow(sigma, 6) / pow(cutoff, 7)) - 2*(pow(sigma, 12) / pow(cutoff, 13))); 
        double potentialEnergyAtCuttoff = 4*epsilon*(pow((sigma/cutoff), 12) - pow((sigma/cutoff), 6));

        double forceNormal = -24*epsilon*((pow(sigma, 6) / pow(distance, 7)) - 2*(pow(sigma, 12) / pow(distance, 13))) - forceCutoff; 
        double forceX = forceNormal * dx / distance;
        double forceY = forceNormal * dy / distance;
        double forceZ = forceNormal * dz / distance;

        particle1.forceX += forceX;
        particle1.forceY += forceY;
        particle1.forceZ += forceZ;
        particle2.forceX -= forceX;
        particle2.forceY -= forceY;
        particle2.forceZ -= forceZ;

        // Stoddard-Ford linear addition
        potentialEnergy += 4*epsilon*(pow((sigma/distance), 12) - pow((sigma/distance), 6)) - potentialEnergyAtCuttoff + forceCutoff * (distance - cutoff);
    }
}

double Collider::GetPotentialEnergy(void){return potentialEnergy;};