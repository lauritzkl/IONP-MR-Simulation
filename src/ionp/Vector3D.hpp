/**
 * @file Vector3D.hpp
 *
 * Legacy two-endpoint vector helper. Norm and dot product use final minus
 * initial position, but several other operations have historical semantics
 * documented below. Simulation diffusion uses Point3D and Particle directly.
 */

/*****************************************************************************
*Class to define a mathematical vector (header .hpp)
*stores initial and final endpoints
*26-8-2021
******************************************************************************/

#ifndef VECTOR3D_HPP_
#define VECTOR3D_HPP_

#include "Point3D.hpp"

class Vector3D{
public:
    Vector3D();
    Vector3D(double x1, double y1, double z1, double x2, double y2, double z2);
    Vector3D(Point3D p1, Point3D p2);

    Vector3D(const Vector3D &v); //copy
    ~Vector3D();

    inline Point3D p_i() const;
    inline Point3D p_f() const;

    void set(double x1, double y1, double z1, double x2, double y2, double z2);
    void print_stdout() const ;

    double norm() const;
    // Legacy atan ratios, not atan2: quadrant information is not recovered.
    double angle_theta() const;
    double angle_phi() const;

    Vector3D &neg();

    double dotProduct(const Vector3D &v) const;

    Vector3D &assign(const Vector3D &v); //self-act operators
    Vector3D &operator=(const Vector3D &v);
    Vector3D &add(const Vector3D &v);
    Vector3D &operator+=(const Vector3D &v);
    Vector3D& subtract(const Vector3D& v);
    Vector3D& operator-=(const Vector3D& v);
    // Scales only the final endpoint; the initial endpoint is unchanged.
    Vector3D &stretch(double a);
    Vector3D &operator*=(double a);

private:
    Point3D m_pi;
    Point3D m_pf;

};

// Legacy comparison behavior: BOTH == and != test left.initial == right.final
// and left.final == right.final. They are not complementary predicates.
bool operator==(const Vector3D &v_left, const Vector3D &v_right); //outside(general) operators
bool operator!=(const Vector3D &v_left, const Vector3D &v_right);
Vector3D operator+(const Vector3D &v_left,const Vector3D &v_right);
Vector3D operator-(const Vector3D& v_left, const Vector3D& v_right);
Vector3D operator*(const Vector3D &v_left, double a);
Vector3D operator*(double a, const Vector3D &v_right);

//inline definitions
inline Vector3D::Vector3D(): m_pi(), m_pf(1.,1.,1.) {} //(0,0,0) (1,1,1)
inline Vector3D::Vector3D(double x1, double y1, double z1, double x2, double y2, double z2): m_pi(x1,y1,z1), m_pf(x2,y2,z2) {}
inline Vector3D::Vector3D(Point3D p1, Point3D p2) : m_pi(p1.x(),p1.y(),p1.z()), m_pf(p2.x(), p2.y(), p2.z()) {}

inline Vector3D::Vector3D(const Vector3D &v): m_pi(v.m_pi), m_pf(v.m_pf){}

inline Point3D Vector3D::p_i() const{
    return m_pi;
}
inline Point3D Vector3D::p_f() const{
    return m_pf;
}

inline Vector3D &Vector3D::operator=(const Vector3D &v){ //return the method result
    return assign(v);
}
inline Vector3D &Vector3D::operator+=(const Vector3D &v){
    return add(v);
}
inline Vector3D& Vector3D::operator-=(const Vector3D& v) {
    return subtract(v);
}
inline Vector3D &Vector3D::operator*=(double a){
    return stretch(a);
}

#endif
