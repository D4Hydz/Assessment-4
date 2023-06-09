#include <iostream>
#include <cmath>
#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

///////////////////////////
// Implement the body class
///////////////////////////
body::body():
pos_(vec(0,0,0)), vel_(vec(0,0,0)), acc_(vec(0,0,0)) 
{
}

body::body(const body &b):
pos_(b.pos_), vel_(b.vel_), acc_(vec(0,0,0))
{
}

body::body(vec pos, vec vel):
pos_(pos), vel_(vel), acc_(vec(0,0,0))
{	
}

void body::set(vec pos, vec vel)
{
	pos_ = pos;
	vel_ = vel;
	acc_ = vec(0,0,0);
}

body &body::operator=(const body &b)
{
	pos_ = b.pos_;
	vel_ = b.vel_;
	acc_ = b.acc_;
	return *this;
}

vec body::get_pos() const
{
	return pos_;
}

vec body::get_vel() const
{
	return vel_;
}

vec body::get_acc() const
{
	return acc_;
}

vec body::get_L() const
{
	return L_;
}

double body::get_ke() const
{
	return ke_;
}

double body::get_gpe() const
{
	return gpe_;
}

void body::set_pos(const vec &v)
{
	pos_ = v;
}

void body::set_vel(const vec &v)
{
	vel_ = v;
}

void body::set_acc(const vec &v)
{
	acc_ = v;
}

void body::set_L(const vec &v)
{
	L_ = v;
}

void body::set_ke(const double E)
{
	ke_ = E;
}

void body::set_gpe(const double E)
{
	gpe_ = E;
}

double body::distance(const body &b) const
{
    vec r = pos_ - b.pos_; // Vector connecting the two bodies
    return r.length(); // Return length of vector
}

vec body::direction(const body &b) const
{
    vec r = pos_ - b.pos_; // Vector connecting the two bodies
    return r.normalize(); // Return normalized vector
}