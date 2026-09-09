/**
 * @file Point3D.hpp
 *
 * Three Cartesian coordinates, used for both positions and displacements.
 * Simulation coordinates use metres. Arithmetic acts component by component;
 * comparisons use exact double equality, with no geometric tolerance.
 */

/*****************************************************************************
*Class to define a 3d point (header .hpp)

*30-8-2021
******************************************************************************/

#ifndef POINT3D_HPP_
#define POINT3D_HPP_

class Point3D{
    public:
    Point3D();
    Point3D(double x, double y, double z);

        Point3D(const Point3D &p); //copy
        ~Point3D() = default;

    inline double x() const;
    inline double y() const;
    inline double z() const;

        inline void set(double x, double y, double z);
    void setX(double x);
    void setY(double x);
    void setZ(double x);
    void print_stdout() const ;

        inline Point3D &assign(const Point3D &p); //self-act operators
    Point3D &operator=(const Point3D &p);
    Point3D &add(const Point3D &p);
    Point3D &operator+=(const Point3D &p);
    Point3D &subtract(const Point3D& p);
    Point3D &operator-=(const Point3D& p);
    Point3D &stretch(double d);
    Point3D &operator*=(double d);

    void placeRandom(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high);

    private:
    double m_x;
    double m_y;
    double m_z;

};

bool operator==(const Point3D &p_left, const Point3D &p_right); //outside(general) operators
bool operator!=(const Point3D &p_left, const Point3D &p_right);
Point3D operator+(const Point3D &p_left,const Point3D &p_right);
Point3D operator-(const Point3D& p_left, const Point3D& p_right);
Point3D operator*(const Point3D &p_left, double d);
Point3D operator*(double d, const Point3D &p_right);

double distance(const Point3D& p_left, const Point3D& p_right); //distance between 2 points

//inline definitions
inline Point3D::Point3D(): m_x(0.), m_y(0.), m_z(0.) {}
inline Point3D::Point3D(double x, double y, double z): m_x(x), m_y(y), m_z(z) {}

inline Point3D::Point3D(const Point3D &p): m_x(p.m_x), m_y(p.m_y), m_z(p.m_z) {}

inline void Point3D::set(double x, double y, double z) {
    m_x = x;
    m_y = y;
    m_z = z;
}

inline Point3D& Point3D::assign(const Point3D& p) {
    if(this != &p) {
        m_x = p.m_x;
        m_y = p.m_y;
        m_z = p.m_z;
    }
    return *this;
}

inline double Point3D::x() const{
    return m_x;
}
inline double Point3D::y() const{
    return m_y;
}
inline double Point3D::z() const{
    return m_z;
}

inline Point3D &Point3D::operator=(const Point3D &p){ //return the method result
    return assign(p);
}
inline Point3D &Point3D::operator+=(const Point3D &p){
    return add(p);
}
inline Point3D& Point3D::operator-=(const Point3D& p) {
    return subtract(p);
}
inline Point3D &Point3D::operator*=(double d){
    return stretch(d);
}

#endif
