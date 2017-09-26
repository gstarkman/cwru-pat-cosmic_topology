#include "chealpix.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex.h>
#include <algorithm>

#include "rotAngles.hpp"
#include "Spherical.hpp"
#include "spline.hpp"

#include "eigenbasis_lin_T3.hpp"
//#include "eigenbasis_sph_T3.hpp"

#ifndef CORRELATION_MATRICES_HPP_INCLUDED
#define CORRELATION_MATRICES_HPP_INCLUDED
class correlation_matrices
{
private:
        int res; //resolution of pixel map
        int lMax; //maximum azimuthal number of the spherical harmonics
        int nMax; //upper limit of the absolute size of the set of integers allowed by each quotient space
public:
        void set_res(int resIn);
        int get_res();
        void set_lMax(int lMaxIn);
        int get_lMax();
        void set_nMax(int nMaxIn);
        int get_nMax();

///     compute the mode-mode correlation matrix of a given fundamental domain
///     inputs:
///         amp: signal amplitude
///         R_LSS: Radius of the Last Scattering Surface
///         k: wave-vector array containing every 3-vector for each set of 3 n's
///         basis: complex-valued array containing the sum of the eigenbasis of the mode functions
///         coeffs: the coefficients of the eigenbasis
///         interpIn: the input needed to perform the interpolation of transfer functions
///         eulAngs: the euler angles of the fundamental domain in the topology frame
        std::vector< std::vector< std::complex<double> > > mode_mode(double amp, double R_LSS, std::vector< std::vector< double > > k,
                                                                        std::vector< std::complex< double > > basis,
                                                                        std::vector< std::vector< double > > interpIn,
                                                                        std::vector< double > eulAngs);

///     transform the mode-mode correlation matrix to the pixel basis
///     input: modeMatIn -> the output of mode_mode or maskCorMat
        std::vector< std::vector< std::complex<double> > > pix_pix(std::vector< std::vector< std::complex<double> > > modeMatIn); //uses res

///     transform one mode-mode correlation matrix to another that is masked
///     inputs:
///         maskIn: array of masking data
///         modeMatIn: the output of mode_mode (technically, could also be the result of another maskCorMat)
        std::vector< std::vector< std::complex<double> > > maskCorMat(std::vector< std::vector< double > > maskIn,
                                                                        std::vector< std::vector< std::complex<double> > > modeMatIn);
};
#endif // CORRELATION_MATRICES_HPP_INCLUDED
