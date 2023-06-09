#ifndef BODY_H
#define BODY_H

#include <iostream>
#include <string>
#include "vector3d.hpp"

// body class definition:
class body
{
public:
	// Initialize with no name, unit mass and zero position, velocity and acceleration vectors
	body();
	// Initialize with another body
    body(const body &b);
    // Initialize explicitly pos, vec
	body(vec pos, vec vel);
	// Set to pos, vec
	void set(vec pos, vec vel);

	// Assign to the values of another body
	body &operator=(const body &b);
	
    // Get the position vector
	vec get_pos() const;

	// Get the velocity vector
	vec get_vel() const;

	// Get the acceleration vector
	vec get_acc() const;

	// Get the angular momentum vector
	vec get_L() const;

	// Get the kinetic energy
	double get_ke() const;

	// Get the gravitational potential energy
	double get_gpe() const;

	// Set the position vector
	void set_name(const std::string n);

	// Set the position vector
	void set_mass(const double m);

	// Set the position vector
	void set_pos(const vec &v);

	// Set the velocity vector
	void set_vel(const vec &v);

	// Set the acceleration vector
	void set_acc(const vec &v);

	// Set the angular momentum vector
	void set_L(const vec &v);

	// Set the kinetic energy
	void set_ke(const double E);

	// Set the gravitational potential energy
	void set_gpe(const double E);

	// Distance between this body and another b
	double distance(const body &b) const;

    // A normalised vector pointing between this body and another b
	vec direction(const body &b) const;



private:
	vec pos_;
	vec vel_;
	vec acc_;
	vec L_;
	double ke_;
	double gpe_;
};

std::ostream &operator<<(std::ostream &os, const body &b); // Overload insertion operator for displaying a body

#endif  // BODY_H