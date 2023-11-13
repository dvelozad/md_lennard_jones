// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
private:
    double mass, radius, x, y, velocityX, velocityY, forceX, forceY;

public:
    void Init(double x0, double y0, double velocityX0, double velocityY0, double mass0, double radius0);
    void Move_r1(double dt);
    void Move_V(double dt);
    void Move_r2(double dt);
    void Draw(void);

    double GetX(void);
    double GetY(void);
    double GetKineticEnergy(void);
    double GetVelocity(void);
    double GetVelocityX(void);

    friend class Collider;
};

#endif
