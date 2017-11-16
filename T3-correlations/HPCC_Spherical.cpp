#include "HPCC_Spherical.hpp"

std::complex<double> cI(0.0, 1.0);

double bigJ(double x, float alpha){
    double J = sqrt(2 / (M_PI * x)) * cos(x - (alpha*M_PI_2) - (M_PI_4));
    if (x == 0.0) return 0.0;
    else return J;
}

double j_l(int l, double x){
    double j = sqrt(M_PI / (2 * x)) * bigJ(x, (l+0.5));
    if (x == 0.0) return 0.0;
    else return j;
}

unsigned int factorial (unsigned int f){
    if (f == 0){
        f = 1;
    }
    else {
        f = f * factorial(f - 1);
    }
    return f;
}

std::complex<double> Y_lm(int l, int m, double polar, double azimuth){
    std::complex<double> ylm = 0.0;
    double plm = 0.0;
    double x = cos(polar);
//    std::cout << " stuff = " << std::sph_legendre(28, 28, 0.0006);
    if (m < 0){
        if (l < 0){ plm = pow(-1,m) * (factorial(abs(l-m)) / factorial(abs(l+m))) * gsl_sf_legendre_sphPlm(abs(l-1), abs(m), x); }
        else { plm = pow(-1,m) * (factorial(abs(l-m)) / factorial(abs(l+m))) * gsl_sf_legendre_sphPlm(l, abs(m), x); }
    }
    else {
        if (l < 0){ plm = gsl_sf_legendre_sphPlm(abs(l-1), abs(m), x); }
        else { plm = gsl_sf_legendre_sphPlm(l, m, x); }
    }
    if (x <= 0){
        plm = std::pow(-1, l+m) * plm;
    }
    ylm = std::pow(-1, m) * plm * std::exp(cI * (m * azimuth));
//    if (std::isnan(std::real(ylm))) std::cout << " plm=" << plm << " azimuth=" << azimuth << " polar=" << polar << std::endl;
    return ylm;
}
