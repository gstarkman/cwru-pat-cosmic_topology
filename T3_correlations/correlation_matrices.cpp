#include "correlation_matrices.hpp"
#include <iostream>
#include <thread>
#include <mutex>
#include <future>
#include <vector>
#include <complex>

std::vector< std::vector< std::complex<double> > > correlation_matrices::mode_mode(double amp, double R_LSS, std::vector< std::vector< double > > k, std::vector< std::complex< double > > basis, std::vector< std::vector< std::complex<double> > > coeffs, std::vector< std::vector< double > > interpIn, std::vector< double > eulAngs){
    std::complex<double> imu(0.0, 1.0);
    int corrsize = (lMax + 1) * (lMax + 1) - 4; //defining the size of the matrix, -4 since we leave out l=0,1
    std::vector< std::vector< std::complex<double> > > modeMat(corrsize, std::vector< std::complex<double> > (corrsize, 0));
    double alpha = eulAngs[0], beta = eulAngs[1]; //gamma = eulAngs[2];
    std::vector<int> lvec(corrsize, 0); // create an array of l's so we can index by nu and nu-prime
    std::vector<int> mvec(corrsize, 0); // create an array of m's ^ditto
    int nu = 0;
    for (int l = 2; l < lMax; l++){
        for (int m = -l; m < l+1; ++m){
            mvec[nu] = m;
            lvec[nu] = l;
            nu++;
        }//m
    }//l
///    initialize interpolation
    std::vector<double> transK = interpIn[0]; // define the x-terms for the cubic spline
    double unitFactr = 1000;
    std::transform(transK.begin(), transK.end(), transK.begin(), std::bind1st(std::multiplies<double>(), unitFactr));
    tk::spline transSpline;
///    correlation matrix
    std::complex<double> corrmatVal = 0.0;
    for (int n = 0; n < k.size(); n++){
        double knorm = sqrt(k[n][0]*k[n][0] + k[n][1]*k[n][1] + k[n][2]*k[n][2]);
    ///    rotation angles
        double theta = rotAngles(1, alpha, k[n]);
        double phi = rotAngles(2, beta, k[n]);
    ///    interpolation
        std::vector<double> transVal;
        std::vector<double> transInterp;
        transInterp.clear();
        for (int ell = 1; ell <= lMax; ell++){
            transVal = interpIn[ell]; // create y-terms for cubic spline; one for each l
            std::transform(transVal.begin(), transVal.end(), transVal.begin(), std::bind1st(std::multiplies<double>(), unitFactr));
            transSpline.set_points(transK, transVal); // initialize the spline
            transInterp.push_back(transSpline(knorm)); // evaluate spline at k
        }//ell
    ///    matrix computation
        for (int v = 0; v < corrsize; v++){
            for (int pv = 0; pv < corrsize; pv++){
                // see eqn
                corrmatVal = std::pow((4 * M_PI), 4) * coeffs[n][lvec[v]] * std::conj(coeffs[n][lvec[pv]]) * std::pow(imu, (lvec[v]-lvec[pv])) *
                                transInterp[lvec[v]] * transInterp[lvec[pv]] *
                                basis[lvec[v]] * std::conj(basis[lvec[pv]]) *
                                j_l(lvec[v], knorm * R_LSS) * std::conj(Y_lm(lvec[v], mvec[v], theta, phi)) *
                                j_l(lvec[pv], knorm * R_LSS) * Y_lm(lvec[pv], mvec[pv], theta, phi);
                modeMat[v][pv] += corrmatVal; //perform sum
            }//pv
        }//v
    }//n
    return modeMat;
}

std::vector< std::vector< std::complex<double> > > correlation_matrices::pix_pix(std::vector< std::vector< std::complex<double> > > modeMatIn){
    int corrsize = (lMax + 1) * (lMax + 1) - 4;
    long nside = pow(2, res); // compute nside of desired map specified by choice of resolution
    long nPix = 12 * nside*nside; // number of pixels in the specified pixel map
    std::vector< std::vector< std::complex<double> > > pixMat(nPix, std::vector< std::complex<double> > (nPix, 0));
    std::vector<int> lvec(corrsize, 0); //ditto, see mode_mode
    std::vector<int> mvec(corrsize, 0);
    int v = 0;
    for (int l = 2; l <= lMax; l++){
        for (int m = -l; m < l+1; ++m){
            mvec[v] = m;
            lvec[v] = l;
            v++;
        }//m
    }//l
    for (int v = 0; v < corrsize; v++){
        for (int pv = 0; pv < corrsize; pv++){
            for (long p1 = 0; p1 < nPix; p1++){
                double phi1, theta1, phi2, theta2; //in healpix, phi is polar and theta is azimuthal so I've switched them so theta is polar in Y_lm and so on
//                pix2ang_nest(nside, p1, &phi1, &theta1);
                pix2ang_ring(nside, p1, &phi1, &theta1);
                for (long p2 = 0; p2 < nPix; p2++){
//                    pix2ang_nest(nside, p2, &phi2, &theta2);
                    pix2ang_ring(nside, p2, &phi2, &theta2);
                    pixMat[p1][p2] += modeMatIn[v][pv] * Y_lm(lvec[v], mvec[v], theta1, phi1) * conj(Y_lm(lvec[pv], mvec[pv], theta2, phi2)); // perform sum
                }//p2 = std::complex<double>(1.0, 1.0);//
            }//p1
        }//pv
    }//v
    return pixMat;
}//uses res, nMax

std::vector< std::vector< std::complex<double> > > correlation_matrices::maskCorMat(std::vector< std::vector<double> > maskingIn, std::vector< std::vector< std::complex<double> > > modeMatIn){
    int corrsize = (lMax + 1) * (lMax + 1) - 4;
    std::vector< std::vector< std::complex<double> > > masking(corrsize, std::vector< std::complex<double> >(corrsize, 0));
    std::vector< std::vector< std::complex<double> > > maskModeMat(corrsize, std::vector< std::complex<double> >(corrsize, 0));
    std::vector<int> lvec(corrsize, 0); //ditto, see mode_mode
    std::vector<int> mvec(corrsize, 0);
    int v = 0;
    for (int l = 2; l <= lMax; l++){
        for (int m = -l; m < l+1; ++m){
            mvec[v] = m;
            lvec[v] = l;
            v++;
        }//m
    }//l
    /// transform the pixel basis masking matrix to the mode basis masking matrix
    for (int p1 = 0; p1 < maskingIn.size(); p1++){
        double phi1, theta1, phi2, theta2; //in healpix, phi is polar and theta is azimuthal so I've switched them so theta is polar in Y_lm and so on
        pix2ang_ring(1024, p1, &phi1, &theta1); // pix2ang_nest is also an option for a nested pixelization scheme here
        for (int p2 = 0; p2 < maskingIn[0].size(); p2++){
            pix2ang_ring(1024, p2, &phi2, &theta2); //ditto
            for (int v = 0; v < corrsize; v++){
                for (int pv = 0; pv < corrsize; pv++){
                    masking[v][pv] += maskingIn[p1][p2] * Y_lm(lvec[v], mvec[v], theta1, phi1) * conj(Y_lm(lvec[pv], mvec[pv], theta2, phi2));
                }//pv
            }//v
        }//p2
    }//p1
    /// transform the mode-mode matrix to the cut sky mode-mode matrix
    for (int v = 0; v < corrsize; v++){
        for (int pv = 0; pv < corrsize; pv++){
            for (int u = 0; u < corrsize; u++){ // u transformed mode index i.e. v -> u; result is also (l,m) x (l',m')
                for (int pu = 0; pu < corrsize; pu++){
                    maskModeMat[u][pu] += modeMatIn[v][pv] * masking[u][v] * std::conj(masking[pu][pv]);
                }//pu
            }//u
        }//pv
    }//v
    return maskModeMat;
}

void correlation_matrices::set_res(int resIn){
    res = resIn;
}

int correlation_matrices::get_res(){
    return res;
}

void correlation_matrices::set_lMax(int lMaxIn){
    lMax = lMaxIn;
}

int correlation_matrices::get_lMax(){
    return lMax;
}

void correlation_matrices::set_nMax(int nMaxIn){
    nMax = nMaxIn;
}

int correlation_matrices::get_nMax(){
    return nMax;
}

