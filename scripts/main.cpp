// main.cpp
#include <iostream>
#include <cmath>
#include "Random64.h"
#include "Particle.h"
#include "Collider.h"
#include "Constants.h"

using namespace std;

void StartAnimation(void) {
    cout << "set terminal gif animate delay 10" << endl; // Set terminal to gif, animate, with a delay
    cout << "set output 'simulation.gif'" << endl; // Output file name
    cout << "unset key" << endl;
    cout << "set xrange [-50:150]" << endl;
    cout << "set yrange [-50:150]" << endl;
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

    StartAnimation();

    collider.Init();

    //PAREDES
    //Pared izquierda
    //particles[N].Init(-10000,Ly/2,0,0,0.01,10000); 
    //Pared derecha
    //particles[N+1].Init(Lx+10000,Ly/2,0,0,0.01,10000); 
    //Pared abajo
    //particles[N+2].Init(Lx/2,-10000,0,0,0.01,10000); 
    //Pared arriba
    //particles[N+3].Init(Lx/2,Ly+10000,0,0,0.01,10000); 


    // Initialize particles
    dx = Lx / (Nx + 1); dy = Ly / (Ny + 1);
    radius = dx / 4;
    //if (dx / 3 > dy / 3) radius = dy / 3;

    for (i = 0; i < N; i++) {
        alpha = 2 * M_PI * randomGenerator.r();
        particles[i].Init(((i % Nx) + 1) * dx, ((i / Ny) + 1) * dy, InitialVelocity * cos(alpha), InitialVelocity * sin(alpha), 10, radius);
    }

    // Main loop
    for (time = drawTime = 0; time < 5000; time += DeltaT, drawTime += DeltaT) {
        // Drawing (optional)

        //if (drawTime > 20 / 120.0) {
        if (int(drawTime) % 10 == 0 && drawTime - int(drawTime) < DeltaT)  {
            StartFrame();
            for (i = 0; i < N; i++) particles[i].Draw();
            EndFrame();
            drawTime = 0;
        }


        // Calculate energies
        kineticEnergy = 0;
        for (i = 0; i < N; i++) kineticEnergy += particles[i].GetKineticEnergy();
        potentialEnergy = collider.GetPotentialEnergy();

        // Uncomment for histogram data
        // if (time > 2)
        //     for (i = 0; i < N; i++) cout << particles[i].GetVelocityX() << endl;

        // Optimized Verlet Velocity
        for (i = 0; i < N; i++) particles[i].Move_r1(DeltaT);
        collider.CalculateForces(particles);
        for (i = 0; i < N; i++) {
            particles[i].Move_V(DeltaT);
            particles[i].Move_r2(DeltaT);
        }
        collider.CalculateForces(particles);
        for (i = 0; i < N; i++) {
            particles[i].Move_V(DeltaT);
            particles[i].Move_r1(DeltaT);
        }
    }

    cout << "set output" << endl; // Close the output file


    return 0;
}
