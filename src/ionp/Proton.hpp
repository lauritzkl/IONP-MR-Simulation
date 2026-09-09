/**
 * @file Proton.hpp
 *
 * Diffusing spin state: Particle position, phase in radians, and cell
 * membership. MRsequence updates this state during diffusion and refocusing.
 * Cellnum is meaningful together with inCell(), not as an independent flag.
 */

/*****************************************************************************
*Class to define a proton (header .hpp)
* derived from particle

*2-9-2021
******************************************************************************/

#ifndef PROTON_HPP_
#define PROTON_HPP_

#include "Point3D.hpp"
#include "Particle.hpp"

class Proton : public Particle {
public:

    inline Proton(double x, double y, double z, double t, double m);
    inline Proton(const Point3D& p, double t, double m);
    inline Proton(double t, double m);

    inline void setPhase(double p) {
        m_phase = p;
    }

    // Historical getter name: the zero-argument overload READS phase in radians.
    inline double setPhase() {
        return m_phase;
    }

    inline void inCell(bool in) {
        m_inCell = in;
    }

    inline bool inCell() const {
        return m_inCell;
    }

    inline void Cellnum(int numCell) {
        m_Cellnum = numCell;
    }

    inline int Cellnum() const {
        return m_Cellnum;
    }

private:
    double m_phase;
    bool m_inCell;
    int m_Cellnum;

};

inline Proton::Proton(double x, double y, double z, double t, double m) : Particle(x,y,z,t), m_phase(m), m_inCell(), m_Cellnum() {}
inline Proton::Proton(const Point3D& p, double t, double m) : Particle(p,t), m_phase(m), m_inCell(), m_Cellnum() {}
inline Proton::Proton(double t, double m) : Particle(0.,0.,0.,t), m_phase(m), m_inCell(), m_Cellnum() {}

#endif
