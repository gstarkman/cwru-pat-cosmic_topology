#include <stdlib.h>
#include <cmath>
#include <complex.h>
#include <vector>
#include <algorithm>

#ifndef ROTANGLES_HPP_INCLUDED
#define ROTANGLES_HPP_INCLUDED

/// compute the angles of rotation needed to transform from the observer frame to the topology frame
///     inputs: angchoice - if 1, return polar rotation angle; if 2, return azimuthal rotation angl
///     euler: the euler angle needed to define the desired rotation
///     vec: the vector we'd like to rotate into the topology frame
double rotAngles(int angchoice, double euler, std::vector<double> vec);

#endif // ROTANGLES_HPP_INCLUDED
