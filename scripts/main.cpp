// main.cpp
#include <iostream>
#include <vector>
#include <array>
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
    std::cout << "Simulation label : " + simulationLabel << endl;
    std::cout << "Box dimension L : " << L << endl;
    std::cout << "Number of particles : " << N << endl;

    Particle particles[N];
    Collider collider;
    Crandom randomGenerator(0);
    double time, drawTime, dx, dy, dz, radius, alpha, kineticEnergy, potentialEnergy, T_current;
    int i;


    if (L / 2 < cutoff) {
        cerr << "Error : The cutoff distance is greater than half the box dimension L / 2 :" << L / 2 << " - Rc : " << cutoff << endl;
        return 1;
    }

    // To save energies
    std::string energyFilename = "../output_files/" + simulationLabel + "_energy_data.txt";
    std::ofstream outFile_energy(energyFilename);
    if (!outFile_energy) {
        cerr << "Error opening file for writing" << endl;
        return 1;
    }

    // To save energies
    std::string temperatureFilename = "../output_files/" + simulationLabel + "_temperature_data.txt";
    std::ofstream outFile_temperature(temperatureFilename);
    if (!outFile_temperature) {
        cerr << "Error opening file for writing" << endl;
        return 1;
    }

    // To save positions
    std::string positionsFilename = "../output_files/" + simulationLabel + "_positions_data.txt";
    std::ofstream outFile_positions(positionsFilename);
    if (!outFile_positions) {
        cerr << "Error opening file for writing" << endl;
        return 1; 
    }

    // To save RDF
    std::string outFile_RDF = "../output_files/" + simulationLabel + "rdf_results.txt";

    // Intit collider
    collider.Init();

    /*    // Initialize particles
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
    }*/

    int unitCellsPerSide = std::cbrt(N / 4);
    double a = Lx / unitCellsPerSide;

    std::vector<std::array<double, 3>> velocities(N);
    std::vector<std::array<double, 3>> positions(N); // Store initial positions

    // Assign random velocities and positions
    double totalVx = 0, totalVy = 0, totalVz = 0;
    int particleIndex = 0;
    for (int ix = 0; ix < unitCellsPerSide; ix++) {
        for (int iy = 0; iy < unitCellsPerSide; iy++) {
            for (int iz = 0; iz < unitCellsPerSide; iz++) {
                std::vector<std::array<double, 3>> unitCellPositions = {
                    {ix * a, iy * a, iz * a},
                    {(ix + 0.5) * a, (iy + 0.5) * a, iz * a},
                    {ix * a, (iy + 0.5) * a, (iz + 0.5) * a},
                    {(ix + 0.5) * a, iy * a, (iz + 0.5) * a}
                };

                for (auto& pos : unitCellPositions) {
                    if (particleIndex < N) {
                        double theta = 2 * M_PI * randomGenerator.r();
                        double phi = acos(2 * randomGenerator.r() - 1);
                        double randomInitialVelocity = randomGenerator.r() * InitialVelocity;

                        double velocityX0 = randomInitialVelocity * sin(phi) * cos(theta);
                        double velocityY0 = randomInitialVelocity * sin(phi) * sin(theta);
                        double velocityZ0 = randomInitialVelocity * cos(phi);

                        velocities[particleIndex] = {velocityX0, velocityY0, velocityZ0};
                        positions[particleIndex] = pos;

                        totalVx += velocityX0;
                        totalVy += velocityY0;
                        totalVz += velocityZ0;

                        particleIndex++;
                    }
                }
            }
        }
    }

    // Adjust velocities to ensure zero average
    double avgVx = totalVx / N;
    double avgVy = totalVy / N;
    double avgVz = totalVz / N;

    for (int i = 0; i < N; i++) {
        particles[i].Init(
            positions[i][0], positions[i][1], positions[i][2],
            velocities[i][0] - avgVx, velocities[i][1] - avgVy, velocities[i][2] - avgVz,
            defaultMass, radius);
    }


    // Steepest Descent Energy Minimization
    for (int step = 0; step < minimizationSteps; step++) {
        collider.CalculateForces(particles); 

        for (i = 0; i < N; i++) {
            particles[i].MinimizeEnergy(minimizationStepSize);
        }
    }

    cout << "Minimization step done" << endl;
    collider.CalculateForces(particles);
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


        double maxDisplacement = 0.0;
        bool shouldUpdate = false;

        for (i = 0; i < N; i++) {
            double displacement = particles[i].CalculateDisplacement();
            maxDisplacement = std::max(maxDisplacement, displacement);

            if (maxDisplacement > (buffer)) {
                shouldUpdate = true;
                break;
            }
        }

        if (shouldUpdate) {
            for (i = 0; i < N; i++) {
                particles[i].UpdateNeighborList(particles);
            }
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

        for (i = 0; i < N; i++) {
            particles[i].Move_r1(dt, 1);
        }

        collider.CalculateForces(particles);
        for (i = 0; i < N; i++) {
            particles[i].Move_V(dt, 0.5);
            //particles[i].Move_r1(dt, 0.5);

            // Update velocities with Langevin thermostat
            particles[i].UpdateVelocity(dt, Gamma, T_desired);
            
        }

        /*        //Omelyan PEFRL
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Zeta);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, (1-2*Lambda)/2);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Xi);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, Lambda);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt,1-2*(Xi+Zeta));
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt, Lambda);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Xi);
        collider.CalculateForces(particles);; for(i = 0; i < N; i++) particles[i].Move_V(dt,( 1-2*Lambda)/2);
        for(i = 0; i < N; i++) particles[i].Move_r1(dt, Zeta);*/


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
