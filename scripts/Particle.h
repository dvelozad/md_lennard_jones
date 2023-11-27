// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
private:
    double mass, radius, x, y, velocityX, velocityY, forceX, forceY;

public:
    void Init(double x0, double y0, double velocityX0, double velocityY0, double mass0, double radius0);
    void Move_r1(double dt, double constant);
    void Move_V(double dt, double constant);

    void Move_V_(double dt);
    void Move_r1_(double dt);
    void Move_r2_(double dt);
    
    void Draw(void);

    double GetMass(void);
    double GetX(void);
    double GetY(void);
    double GetKineticEnergy(void);
    double GetVelocity(void);
    double GetVelocityX(void);
    double GetVelocityY(void);

    void ApplyPeriodicBoundaryConditions(double Lx, double Ly);
    void MinimizeEnergy(double dt);

    void RescaleVelocity(double lambda);

    friend class Collider;
};

#endif