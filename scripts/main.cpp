// main.cpp
#include <iostream>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "Random64.h"
#include "Particle.h"
#include "Collider.h"
#include "Constants.h"

using namespace std;


double CalculateCurrentTemperature(Particle *particles, int numParticles) {
    double totalKineticEnergy = 0.0;

    for (int i = 0; i < numParticles; ++i) {
        double velocitySquared = particles[i].GetVelocity() * particles[i].GetVelocity();
        totalKineticEnergy += 0.5 * particles[i].GetMass() * velocitySquared;
    }

    double temperature = (2.0 / (3.0 * numParticles * kB)) * totalKineticEnergy;
    return temperature;
}

// Main function
int main() {
    Particle particles[N];
    Collider collider;
    Crandom randomGenerator(0);
    double time, drawTime, dx, dy, radius, alpha, kineticEnergy, potentialEnergy;  
    int i;

    // To save energies
    std::ofstream outFile("../output_files/energy_data.txt");
    if (!outFile) {
        cerr << "Error opening file for writing" << endl;
        return 1; // or handle the error as you see fit
    }

    // To save energies
    std::ofstream outFile_positions("../output_files/positions_data.txt");
    if (!outFile_positions) {
        cerr << "Error opening file for writing" << endl;
        return 1; // or handle the error as you see fit
    }

    double defaultMass = 1;
    double totalTime = 10000;
    int timeFrame = 5;

    //StartAnimation();

    collider.Init();

    // Initialize particles
    dx = Lx / (Nx + 1); dy = Ly / (Ny + 1);
    radius = dx / 8;
    //if (dx / 3 > dy / 3) radius = dy / 3;

    for (i = 0; i < N; i++) {
        alpha = 2 * M_PI * randomGenerator.r();
        double randomInitialVelocity =  InitialVelocity * randomGenerator.r();
        particles[i].Init(((i % Nx) + 1) * dx, ((i / Ny) + 1) * dy, randomInitialVelocity * cos(alpha), randomInitialVelocity * sin(alpha), defaultMass, radius);
        //particles[i].Init(((i % Nx) + 1 + 3*randomGenerator.r()) * dx, ((i / Ny) + 1 + randomGenerator.r()) * dy , randomInitialVelocity * cos(alpha), randomInitialVelocity * sin(alpha), defaultMass, radius);
    }

    // Steepest Descent Energy Minimization
    double minimizationStepSize = 0.0001;  
    int minimizationSteps = 100000;        

    for (int step = 0; step < minimizationSteps; step++) {
        collider.CalculateForces(particles); 

        for (i = 0; i < N; i++) {
            particles[i].MinimizeEnergy(minimizationStepSize);
        }
    }

    cout << "Minimization step done" << endl;

    for (time = drawTime = 0; time < totalTime; time += dt, drawTime += dt) {
        // Drawing (optional)

        //if (drawTime > 20 / 120.0) {
        if (int(drawTime) % timeFrame == 0 && drawTime - int(drawTime) < dt)  {

            for (i = 0; i < N; i++){
                outFile_positions << i << " " << time << " " << particles[i].GetX() << " " << particles[i].GetY() << endl;
            }
        }

        // Uncomment for histogram data
        // if (time > 2)
        //     for (i = 0; i < N; i++) cout << particles[i].GetVelocityX() << endl;

        /*        // Optimized Verlet Velocity - Old algorithm
        for (i = 0; i < N; i++) particles[i].Move_r1_(dt);

        collider.CalculateForces(particles);

        for (i = 0; i < N; i++) {
            particles[i].Move_V_(dt);
            particles[i].Move_r2_(dt);
        }
        
        collider.CalculateForces(particles);

        for (i = 0; i < N; i++) {
            particles[i].Move_V_(dt);
            particles[i].Move_r1_(dt);
        }*/


        // Calculate energies
        if (int(drawTime) % timeFrame == 0 && drawTime - int(drawTime) < dt)  {
            kineticEnergy = 0;
            for (i = 0; i < N; i++) kineticEnergy += particles[i].GetKineticEnergy();
            potentialEnergy = collider.GetPotentialEnergy();
            outFile << time << " " << kineticEnergy + potentialEnergy << endl;
        }

        //Omelyan PEFRL
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Zeta);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, (1-2*Lambda)/2);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Xi);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, Lambda);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt,1-2*(Xi+Zeta));
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, Lambda);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Xi);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt,( 1-2*Lambda)/2);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Zeta);


        // After calculating the current temperature of the system
        double T_current = CalculateCurrentTemperature(particles, N);
        double lambda = sqrt(T_desired / T_current);

        if(T_current == 0){
            cout << "Current temperature is zero!!" << endl;
        }

        // Rescale velocities
        for (int i = 0; i < N; i++) {
            particles[i].RescaleVelocity(lambda);
        }
    }

    return 0;
}
