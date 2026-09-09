#include <cmath>
#include <cassert>
#include <iostream>

#include "Point3D.hpp"
#include "Proton.hpp"
#include "MRsequence.hpp"

// Forward declaration of the helper defined in MRsequence.cpp
bool reflect_from_sphere(const Point3D& sphere_center, double sphere_rad, Proton& proton, Point3D& last_pos, double current_time, std::vector<TrajectoryEvent>* event_log);

// Simple deterministic check: a proton moving along +x toward the origin should
// hit the sphere at x = -r and reflect back along -x, ending at x = start - 2*r
// minus the remaining path length. This validates specular reflection and the
// last_pos update.
int main() {
	// Sphere centered at origin with radius 1
	const Point3D center(0.0, 0.0, 0.0);
	const double radius = 1.0;

	// Proton travels from x=-2 to x=+2 (length 4), hits at x=-1, remaining length 3 -> final x = -1 - 3 = -4
	Point3D last_pos(-2.0, 0.0, 0.0);
	Proton proton(0.0, 0.0); // phase/time ctor
	proton.setXt(Point3D(2.0, 0.0, 0.0));

	const bool hit = reflect_from_sphere(center, radius, proton, last_pos, 0.0, nullptr);
	if (!hit) {
		std::cerr << "Expected an intersection with the sphere\n";
		return 1;
	}

	// last_pos is nudged just outside the hit point to avoid repeated self-intersections.
	assert(std::abs(last_pos.x() + 1.0) < 2e-6);
	assert(std::abs(last_pos.y()) < 1e-12);
	assert(std::abs(last_pos.z()) < 1e-12);

	// Reflected endpoint should be at x=-4 after continuing remaining distance in reflected dir
	assert(std::abs(proton.xt().x() + 4.0) < 1e-12);
	assert(std::abs(proton.xt().y()) < 1e-12);
	assert(std::abs(proton.xt().z()) < 1e-12);

	std::cout << "Elastic reflection test passed\n";
	return 0;
}
