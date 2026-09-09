/**
 * @file Aggregate.hpp
 *
 * Aggregate geometry: a spherical placement region containing a requested
 * number of IONPs. Positions and radii use metres. Particle supplies the
 * reference/current positions; SimSpace performs the actual IONP placement.
 */

/*****************************************************************************
*Class to define an aggregate (header .hpp)

*2-9-2021
******************************************************************************/

#ifndef AGGREGATE_HPP_
#define AGGREGATE_HPP_

#include "Point3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"

class Aggregate : public Particle{
public:
    Aggregate(double x = 0., double y = 0., double z = 0., double t=0.,  double r = 0, const int n = 1);
    Aggregate(const Point3D& p, double t = 0., double r = 0., int n = 1);
    Aggregate(double r = 0., int n = 1);
    Aggregate(int n = 1);

    inline double diffTime() const;
    inline double radius() const;
    inline int nparticles() const;

private:
    inline void setDiffTime(double t); // Hides the base setter; this class retains separate time metadata.
    void diffuse(double d, double dt);

    double m_diff_time = 0;
    double m_radious;
    int m_nparticles;
};

//inline definitions
inline Aggregate::Aggregate(double x, double y, double z, double t, double r, int n) : Particle(x,y,z,t), m_radious(r), m_nparticles(n)  {}
inline Aggregate::Aggregate(const Point3D& p, double t, double r, int n) : Particle(p, t), m_radious(r), m_nparticles(n) {}
inline Aggregate::Aggregate(double r, int n) : Particle(0.,0.,0.,0.), m_radious(r), m_nparticles(n) {}
inline Aggregate::Aggregate(int n) : Particle(0.,0.,0.,0.), m_radious(), m_nparticles(n) {}

inline double Aggregate::diffTime() const {
    return m_diff_time;
}
inline double Aggregate::radius() const {
    return m_radious;
}
inline int Aggregate::nparticles() const {
    return m_nparticles;
}

inline void Aggregate::setDiffTime(double t) {
    m_diff_time = t;
}

#endif
