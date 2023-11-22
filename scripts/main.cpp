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

void StartAnimation(void) {
    cout << "set terminal gif animate delay 10" << endl; // Set terminal to gif, animate, with a delay
    cout << "set output '../output_files/simulation.gif'" << endl; // Output file name
    cout << "unset key" << endl;
    cout << "set xrange [-5:55]" << endl;
    cout << "set yrange [-5:55]" << endl;
    cout << "set size ratio -1" << endl;
    cout << "set parametric" << endl;
    cout << "set trange [0:7]" << endl;
    cout << "set isosamples 12" << endl;
}


// Animation functions continued
void StartFrame(void) {
    cout << "plot 0,0 ";
    // Add more plot commands as needed
}

void EndFrame(void) {
    cout << endl;
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

    double totalTime = 100000;

    StartAnimation();

    collider.Init();

    // Initialize particles
    dx = Lx / (Nx + 1); dy = Ly / (Ny + 1);
    radius = dx / 8;
    //if (dx / 3 > dy / 3) radius = dy / 3;

    for (i = 0; i < N; i++) {
        alpha = 2 * M_PI * randomGenerator.r();
        particles[i].Init(((i % Nx) + 1) * dx, ((i / Ny) + 1) * dy, InitialVelocity * cos(alpha), InitialVelocity * sin(alpha), 10, radius);
    }

    // Main loop
    for (time = drawTime = 0; time < totalTime; time += dt, drawTime += dt) {
        // Drawing (optional)

        //if (drawTime > 20 / 120.0) {
        if (int(drawTime) % 70 == 0 && drawTime - int(drawTime) < dt)  {

            for (i = 0; i < N; i++){
                outFile_positions << i << " " << time << " " << particles[i].GetX() << " " << particles[i].GetY() << endl;
            }
            /*            StartFrame();
            for (i = 0; i < N; i++) particles[i].Draw();
            EndFrame();
            drawTime = 0;*/
        }

        // Calculate energies

        if (int(drawTime) % 70 == 0 && drawTime - int(drawTime) < dt)  {
            kineticEnergy = 0;
            for (i = 0; i < N; i++) kineticEnergy += particles[i].GetKineticEnergy();
            potentialEnergy = collider.GetPotentialEnergy();
            outFile << time << " " << kineticEnergy + potentialEnergy << endl;
        }

        // Uncomment for histogram data
        // if (time > 2)
        //     for (i = 0; i < N; i++) cout << particles[i].GetVelocityX() << endl;

        /*
        // Optimized Verlet Velocity - Old algorithm
        for (i = 0; i < N; i++) particles[i].Move_r1(dt);

        collider.CalculateForces(particles);

        for (i = 0; i < N; i++) {
            particles[i].Move_V(dt);
            particles[i].Move_r2(dt);
        }
        
        collider.CalculateForces(particles);

        for (i = 0; i < N; i++) {
            particles[i].Move_V(dt);
            particles[i].Move_r1(dt);
        }
        */

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
    }

    cout << "set output" << endl; // Close the output file


    return 0;
}
