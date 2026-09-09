/**
 * @file Point3D.cpp
 *
 * Component-wise Point3D operations. Assignment and compound operators mutate
 * the receiver; free arithmetic operators return a modified copy.
 */

/*****************************************************************************
*Class to define a 3d point (cpp)

*30-8-2021
******************************************************************************/

#include <iostream>
#include <cmath>
#include "mptmacros.h"
#include "Point3D.hpp"
#include "RandNumberBetween.hpp"

void Point3D::setX(double x){
    m_x = x;
}

void Point3D::setY(double x){
    m_y = x;
}

void Point3D::setZ(double x){
    m_z = x;
}

void Point3D::print_stdout() const{
    std::cout<<"Coordinates (x,y,z) = (" <<m_x<<","<<m_y<<","<<m_z<<")"<<std::endl;
}

Point3D &Point3D::add(const Point3D &p){
    m_x += p.x();
    m_y += p.y();
    m_z += p.z();

    return *this;
}

Point3D& Point3D::subtract(const Point3D& p) {
    m_x -= p.x();
    m_y -= p.y();
    m_z -= p.z();

    return *this;
}

Point3D &Point3D::stretch(double d){
    m_x *= d;
    m_y *= d;
    m_z *= d;

    return *this;
}

void Point3D::placeRandom(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high) { //different limits

    RandNumberBetween x_rand_numb(x_low, x_high);
    RandNumberBetween y_rand_numb(y_low, y_high);
    RandNumberBetween z_rand_numb(z_low, z_high);

  m_x = x_rand_numb();
  m_y = y_rand_numb();
  m_z = z_rand_numb();
}

bool operator==(const Point3D &p_left,const Point3D &p_right){
    return p_left.x()==p_right.x() && p_left.y()==p_right.y() && p_left.z()==p_right.z();
}

bool operator!=(const Point3D &p_left,const Point3D &p_right){
    return p_left.x()!=p_right.x() || p_left.y()!=p_right.y() || p_left.z()!=p_right.z();
}

Point3D operator+(const Point3D &p_left,const Point3D &p_right){ //left to right addition by using class methods
    return Point3D(p_left).add(p_right);
}

Point3D operator-(const Point3D& p_left, const Point3D& p_right) { //left to right addition by using class methods
    return Point3D(p_left).subtract(p_right);
}

Point3D operator*(const Point3D &p_left, double d){
    return Point3D(p_left).stretch(d);
}

Point3D operator*(double d, const Point3D &p_right){
    return Point3D(p_right).stretch(d);
}

double distance(const Point3D& p_left, const Point3D& p_right) { //distance between 2 points
    double dist = std::sqrt(square(p_left.x() - p_right.x()) + square(p_left.y() - p_right.y()) + square(p_left.z() - p_right.z()) );

    return dist;
}
