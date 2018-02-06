#ifndef HPCC_SPHERICAL_HPP_INCLUDED
#define HPCC_SPHERICAL_HPP_INCLUDED

#include <complex>
#include <stdlib.h>
#include <vector>
#include <gsl_sf_result.h>
#include <gsl_sf_legendre.h>

/// these are the typical mathematical functions for expansion in spherical harmonics
double j_l(int l, double x);
std::complex<double> Y_lm(int l, int m, double polar, double azimuth);

#endif // HPCC_SPHERICAL_HPP_INCLUDED
