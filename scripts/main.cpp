// main.cpp
#include <iostream>
#include <vector>
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


void ComputeAndWriteRDF(Particle *particles, double maxDistance, int numBins, const std::string& filename) {
    std::vector<int> bins(numBins, 0);
    double binSize = maxDistance / numBins;
    std::vector<double> rdf(numBins, 0.0);

    // Calculate RDF
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double distance = particles[i].DistanceTo(particles[j]);

            int binIndex = distance / binSize;
            if (binIndex < numBins) {
                bins[binIndex]++;
            }
        }
    }

    // write
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Normalize RDF 
    double density = N / (Lx * Ly); 
    for (int i = 0; i < numBins; i++) {
        double r1 = i * binSize;
        double r2 = r1 + binSize;

        // 2D case shell volume
        //double shellVolume = M_PI * (r2 * r2 - r1 * r1);  

        // 3D case shell volume
        double shellVolume = (4.0 / 3.0) * M_PI * (r2* r2 * r2 - r1 * r1 * r1);  

        rdf[i] = (bins[i] / N) / (density * shellVolume); // Avg number of pairs in that bin div the number of ideal gas partciles at that density

         // Writing the bin center and RDF value
        outFile << (r1 + r2) / 2 << "\t" << rdf[i] << "\n"; 
    }

    outFile.close();
}


// Main function
int main() {
    Particle particles[N];
    Collider collider;
    Crandom randomGenerator(0);
    double time, drawTime, dx, dy, dz, radius, alpha, kineticEnergy, potentialEnergy, T_current;
    int i;

    // To save energies
    std::ofstream outFile_energy("../output_files/energy_data.txt");
    if (!outFile_energy) {
        cerr << "Error opening file for writing" << endl;
        return 1;
    }

    // To save energies
    std::ofstream outFile_temperature("../output_files/temperature_data.txt");
    if (!outFile_temperature) {
        cerr << "Error opening file for writing" << endl;
        return 1;
    }

    // To save positions
    std::ofstream outFile_positions("../output_files/positions_data.txt");
    if (!outFile_positions) {
        cerr << "Error opening file for writing" << endl;
        return 1; 
    }

    // To save RDF
    std::string outFile_RDF = "../output_files/rdf_results.txt";

    // Intit collider
    collider.Init();

    // Initialize particles
    dx = Lx / (Nx + 1); dy = Ly / (Ny + 1); dz = Lz / (Nz + 1);
    radius = dx / 8;
    for (i = 0; i < N; i++) {
        // Generating random angles for spherical symmetry
        double theta = 2 * M_PI * randomGenerator.r(); 
        double phi = acos(2 * randomGenerator.r() - 1); 

        // Random velocity magnitude
        double randomInitialVelocity = InitialVelocity * randomGenerator.r();

        // Velocity components
        double velocityX0 = randomInitialVelocity * sin(phi) * cos(theta);
        double velocityY0 = randomInitialVelocity * sin(phi) * sin(theta);
        double velocityZ0 = randomInitialVelocity * cos(phi);

        // Initialize particle
        particles[i].Init(
            ((i % Nx) + 1) * dx,
            ((i / Nx) % Ny + 1) * dy,
            ((i / (Nx * Ny)) + 1) * dz,
            velocityX0, velocityY0, velocityZ0, 
            defaultMass, radius
        );
        //particles[i].Init(((i % Nx) + 1 + 3*randomGenerator.r()) * dx, ((i / Ny) + 1 + randomGenerator.r()) * dy , randomInitialVelocity * cos(alpha), randomInitialVelocity * sin(alpha), defaultMass, radius);
    }

    // Steepest Descent Energy Minimization
    for (int step = 0; step < minimizationSteps; step++) {
        collider.CalculateForces(particles); 

        for (i = 0; i < N; i++) {
            particles[i].MinimizeEnergy(minimizationStepSize);
        }
    }

    cout << "Minimization step done" << endl;

    for (time = drawTime = 0; time < totalTime; time += dt, drawTime += dt) {

        // Write info
        if (int(drawTime) % timeFrame == 0 && drawTime - int(drawTime) < dt)  {
            // Get positions
            for (i = 0; i < N; i++){
                outFile_positions << i << " " << time << " " << particles[i].GetX() << " " << particles[i].GetY() << " " << particles[i].GetZ() << endl;
            }

            // Get energy
            kineticEnergy = 0;
            for (i = 0; i < N; i++) kineticEnergy += particles[i].GetKineticEnergy();
            potentialEnergy = collider.GetPotentialEnergy();

            // Get temperature
            T_current = CalculateCurrentTemperature(particles, N);

            // Write
            outFile_energy << time << " " << kineticEnergy + potentialEnergy << endl;
            outFile_temperature << time << " " << T_current << endl;
        }       

        // Histogram data
        // if (time > 2)
        //     for (i = 0; i < N; i++) cout << particles[i].GetVelocityX() << endl;

        /*        
        // Optimized Verlet Velocity - Old algorithm
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
    
        // Update velocities with Langevin thermostat
        for (int i = 0; i < N; i++) {
            particles[i].UpdateVelocity(dt, Gamma, T_desired);
        }

        // After calculating the current temperature of the system
        //double T_current = CalculateCurrentTemperature(particles, N);
        //double lambda = sqrt(T_desired / T_current);

        //if(T_current == 0){
        //    cout << "Current temperature is zero!!" << endl;
        //}

        // Rescale velocities
        //for (int i = 0; i < N; i++) {
        //    particles[i].RescaleVelocity(lambda);
        //}

        //double T_current = CalculateCurrentTemperature(particles, N);

        /*        // Update xi and xidot
        kineticEnergy = 0;
        for (i = 0; i < N; i++) kineticEnergy += particles[i].GetKineticEnergy();

        for (int i = 0; i < N; i++) {
            particles[i].UpdateThermostat(kineticEnergy);
        }*/
    }

    // Write RDF
    //ComputeAndWriteRDF(particles, maxDistance, numBins, outFile_RDF);

    return 0;
}
