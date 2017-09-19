#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <complex>
#include <stdlib.h>
#include <vector>

#ifndef SPHERICAL_HPP_INCLUDED
#define SPHERICAL_HPP_INCLUDED

/// these are the typical mathematical functions for expansion in spherical harmonics
double j_l(int l, double x);
std::complex<double> Y_lm(int l, int m, double polar, double azimuth);

#endif // SPHERICAL_HPP_INCLUDED
