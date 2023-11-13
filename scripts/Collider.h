// Collider.h
#ifndef COLLIDER_H
#define COLLIDER_H

#include "Particle.h"

class Collider {
private:
    double potentialEnergy;

public:
    void Init(void);
    void CalculateForces(Particle *particles);
    void Collide(Particle &particle1, Particle &particle2);
    double GetPotentialEnergy(void);
};

#endif
