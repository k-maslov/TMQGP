
#include <Particle.h>

Particle::Particle(double m, double * erange, int dimE, double * qrange, int dimQ, bool stat, double eps, double d)
{
    this->dimE = dimE;
    this->dimQ = dimQ;

    this->erange = erange;
    this->qrange = qrange;
    
    
}
