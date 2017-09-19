#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>

#include "Spherical.hpp"

#ifndef EIGENBASIS_LIN_T3_HPP_INCLUDED
#define EIGENBASIS_LIN_T3_HPP_INCLUDED

class eigenbasis_lin_T3
{
private:
        int nMax; //upper limit of the absolute size of the set of integers allowed by each quotient space
        int domIndx; //index of the quotient space the computation is to be performed for
        std::vector<double> angles = std::vector<double>(3); //angles parameterizing the fundamental domain
        std::vector<double> lengths = std::vector<double>(3); //lengths parameterizing the fundamental domain
public:
/// 2D arrays listing all necessary integers for a chosen domIndx
        struct retNs {
            std::vector< std::vector<int> > ns; //1st dimension is 3; 2nd dimension depends on choice of nMax
            std::vector< std::vector<int> > nsTr; //1st dim depends on choice of nMax; 2nd dimension is 3
        };
/// storage for all quantities used to compute the eigenmode basis that are also needed for CorrMat
        struct retVals {
            std::vector< std::vector<double> > k;           //contains the wave-vectors. size depends on nMax
            std::vector< std::complex<double> > basis;      //contains the eigenmode basis sum. size depends on nMax
            std::vector< std::vector< std::complex<double> > > coeffs;      //contains basis coefficients
        };

        void set_nMax(int nMaxIn);
        int get_nMax();
        void set_domIndx(int domIndxIn);
        int get_domIndx();
        void set_angles(std::vector<double> anglesIn);
        std::vector<double> get_angles();
        void set_lengths(std::vector<double> lengthsIn);
        std::vector<double> get_lengths();

///     compute each integer allowed by the quotient space up to nMax and store both the set of three sets of allowed n's
///     and a list of n-vectors, which are 3-vectors
        retNs domainNvec();

///     compute the terms needed to perform the eigenmode basis summation and store them as well as the summation result
        retVals eigBasis(double Rlss, std::vector<double> rT);
};


#endif // EIGENBASIS_LIN_T3_HPP_INCLUDED
