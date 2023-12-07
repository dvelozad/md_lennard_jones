// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle {
private:
    double mass, radius, x, y, z, velocityX, velocityY, velocityZ, forceX, forceY, forceZ;
    double xi = 0.0; // Thermostat variable
    double xidot = 0.0; // Time derivative of the thermostat variable

    double forceOX = 0; double forceOY = 0; double forceOZ = 0;

    double lastX, lastY, lastZ; // Store the positions at last neighbor list update

public:
    void Init(double x0, double y0, double z0, double velocityX0, double velocityY0, double velocityZ0, double mass0, double radius0);
    void Move_r1(double dt, double constant);
    void Move_V(double dt, double constant);

    void Move_V_(double dt);
    void Move_r1_(double dt);
    void Move_r2_(double dt);
    
    void Draw(void);

    double GetMass(void);
    double GetX(void);
    double GetY(void);
    double GetZ(void);
    double GetKineticEnergy(void);
    double GetVelocity(void);
    double GetVelocityX(void);
    double GetVelocityY(void);
    double GetVelocityZ(void);
    double GetNeighborList(void);

    void ApplyPeriodicBoundaryConditions(double Lx, double Ly, double Lz);
    void MinimizeEnergy(double dt);

    void RescaleVelocity(double lambda);

    void UpdateVelocity(double dt, double gamma, double temperature);

    double DistanceTo(Particle &other);

    std::vector<int> neighborList;
    void UpdateNeighborList(Particle *particles);
    double CalculateDisplacement();


    friend class Collider;
};

#endif