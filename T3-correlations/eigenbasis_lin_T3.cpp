#include "eigenbasis_lin_T3.hpp"
#include <iostream>

std::complex<double> imu(0.0, 1.0); //imaginary unit

/**==============================__Wave-Vectors__===================================**/
/// some notation:
/// KEI -> Ith fundamental domain's wave-vector
/// KEI_NJ -> N = J solution

/// inputs:
/// nvec: 3-vector of elements of each allowed set of integers; nvec[0],[1],[2] <-> nx,ny,nz
/// angles: array of 3 tilt angles defining fundamental domain; angles[0],[1],[2] <-> theta, phi ,psi
/// lengths: array of 3 lengths defining fundamental domain; lengths[0],[1],[2] <-> lx,ly,lz -or- l,l,lz

/// outputs:
/// k: 3-vector of 3 wave numbers in each principle direction
std::vector<double> K(std::vector<int> nvec, std::vector<double> angles, std::vector<double> lengths){
    std::vector<double> k(3, 0.0);
    k[0] = 2 * M_PI * ( (nvec[0] / lengths[0]) + (nvec[1] / lengths[1]) * std::cos(angles[0]) / std::sin(angles[0]) );
    k[1] = 2 * M_PI * ( (nvec[1] / lengths[1]) + (nvec[2] / lengths[2]) * std::cos(angles[1]) / std::sin(angles[1]) );
    k[2] = 2 * M_PI * ( (nvec[2] / lengths[2]) + ( (nvec[0] / lengths[0]) + (nvec[1] / lengths[1]) * std::cos(angles[0]) / std::sin(angles[0]) ) * std::cos(angles[2]) / std::sin(angles[2]) );
    return k;
}

/**================================__Coefficients__=================================**/
/// some notation:
/// CEI -> Ith fundamental domain's eigenmode's basis coefficients

/// inputs:
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI
/// nvec: 3-vector of elements of each allowed set of integers

/// outputs:
/// coeffVec: complex-valued 3-vector containing each a_q needed for representing
///           the eigenmodes as linear combinations of the covering space eigenmodes of domain I

std::vector< std::complex<double> > expMatArg(std::vector<double> vec1,
                                              std::vector<double> vec2,
                                              std::vector< std::vector< std::vector<double> > > Mlist){
    std::vector< std::complex<double> > expVec(Mlist.size(), 0);
    for (int i = 0; i < Mlist.size(); i++){
        std::complex<double> expTerm = 1.0;
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                expTerm *= std::exp(imu*vec1[j]*Mlist[i][j][k]*vec2[k]);
            }
        }
        expVec[i] = expTerm;
    }
    return expVec;
}

std::complex<double> sum_EI(std::vector< std::complex<double> > xiterm,
                            std::vector< std::complex<double> > Cf) {
    std::complex<double> esum = 0;
    for (int i = 0; i < Cf.size(); i++){
        esum += Cf[i]*xiterm[i];
    }//i
    return esum;
}

std::vector< std::vector<double> > matPow(std::vector< std::vector<double> > matIn, int pow){ //compute m^n like m.m.m. ... .m for 3x3 m
    if (pow == 0){
        std::vector< std::vector<double> > matOut(3, std::vector<double>(3, 0));
        matOut[0][0] = 1.0; matOut[1][1] = 1.0; matOut[2][2] = 1.0;
        return matOut;
    }
    else if (pow == 1){
        return matIn;
    }
    else{
        std::vector< std::vector<double> > matOut(3, std::vector<double>(3, 0));
        matOut = matIn;
        for (int i = 0; i < 3; i++){ //rows
            for (int j = 0; j < 3; j++){ //columns
                int p = 2;
                while (p < pow){ //repeat for however many powers of matIn you want
                        matOut[i][j] += matOut[i][j] * matIn[i][j];
                        p++;
                }//p
            }//j
        }//i
        return matOut;
    }
}

/**==============================__Normalization_Constants__===================================**/
/// inputs:
/// lengths: 3 vector of fundamental domain side lengths
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI; must match domIndx to NEI as well
/// coeffs: coefficients of the eigenmode's basis; outputs of CEI functions
/// transMats: vector of transformation matrices or products of transfomation matrices

/// output:
/// 1.0 / sqrt(N): The normalization constant for the linear representation of eigenmode functions for the domain I

std::complex<double> expComp (double x, double k, double m1, double m2) {
    return imu * std::exp(imu*k*m1*x) * std::exp(-imu*k*m2*x) / (k * (m1-m2));
};

std::complex<double> NConst_rprism (std::vector<double> lengths,
                                    std::vector<double> kvec,
                                    std::vector< std::complex<double> > Cf,
                                    std::vector< std::vector< std::vector<double> > > Mlist){
    int q_max = Mlist.size();
    std::complex<double> N = 0.0;
    std::complex<double> prod = 1.0;
    for (int q; q < q_max; q++) {
        for (int qp; qp < q_max; qp++) {
            for (int r; r < 3; r++) {
                for (int s; s < 3; s++) {
                    prod *= expComp(lengths[s], kvec[r], Mlist[q][r][s], Mlist[qp][r][s]) -
                            expComp(0, kvec[r], Mlist[q][r][s], Mlist[qp][r][s]);
                }//s
            }//r
            N += Cf[q]* std::conj(Cf[qp]) * prod;
        }//qp
    }//q
    return N;
}

std::complex<double> expComp_1_2 (double x, double y, double k, double m1, double m1p, double m2, double m2p, double m3, double m3p, double b1, double b2) {
    return (-imu * std::exp(imu*k*x*(m1 - m1p + 4 * (m2 - m2p)) + y * (m3 - m3p)) / (k * (m2 - m2p) * (m3 - m3p))) *
           (
                std::exp(imu*k*(b1*(m2 - m2p) + x*(m3 - m3p))) / (k * (m3 - m3p + 4*(m2 - m2p))) +
                std::exp(imu*k*(b2*(m2 - m2p) - x*(m3 - m3p))) / (k * (m3 + m3p - 4*(m2 - m2p)))
           );
};

std::complex<double> NConst_hex (std::vector<double> lengths,
                                 std::vector<double> kvec,
                                 std::vector< std::complex<double> > Cf,
                                 std::vector< std::vector< std::vector<double> > > Mlist){
    int q_max = Mlist.size();
    std::complex<double> N = 0.0;
    double b_AB = lengths[1] / 2;
    double b_BC = b_AB;
    double b_DE = 0.5*(lengths[1] - (lengths[1]*lengths[0]/(0.375*lengths[1] - lengths[0])));
    double b_EF = 3*lengths[1] / (16*lengths[0] - 6*lengths[1]);
    std::vector<double> ulims(3, 0);
    ulims.push_back(0.375*lengths[1]); ulims.push_back(lengths[1]); ulims.push_back(lengths[2]);
    std::vector<double> llims(3, 0);
    llims.push_back(0.25*lengths[1]); llims.push_back(0.0); llims.push_back(0.0);
    std::complex<double> N1qq, N2qq, N3qq;
    for (int q; q < q_max; q++) {
        for (int qp; qp < q_max; qp++) {
            N1qq = 1.0; N2qq = 1.0; N3qq = 1.0;
            for (int r; r < 3; r++) {
                for (int s; s < 3; s++) {
                    N1qq *= expComp_1_2(0.25*lengths[1], lengths[2],
                                        kvec[r],
                                        Mlist[q][r][0], Mlist[qp][r][0],
                                        Mlist[q][r][1], Mlist[qp][r][1],
                                        Mlist[q][r][2], Mlist[qp][r][2],
                                        b_AB, b_BC) -
                            expComp_1_2(0.0, 0.0,
                                        kvec[r],
                                        Mlist[q][r][0], Mlist[qp][r][0],
                                        Mlist[q][r][1], Mlist[qp][r][1],
                                        Mlist[q][r][2], Mlist[qp][r][2],
                                        b_AB, b_BC);
                    N2qq *= expComp_1_2(lengths[0], lengths[2],
                                        kvec[r],
                                        Mlist[q][r][0], Mlist[qp][r][0],
                                        Mlist[q][r][1], Mlist[qp][r][1],
                                        Mlist[q][r][2], Mlist[qp][r][2],
                                        b_EF, b_DE) -
                            expComp_1_2(0.375*lengths[1], 0.0,
                                        kvec[r],
                                        Mlist[q][r][0], Mlist[qp][r][0],
                                        Mlist[q][r][1], Mlist[qp][r][1],
                                        Mlist[q][r][2], Mlist[qp][r][2],
                                        b_EF, b_DE);
                    N3qq *= expComp(ulims[s], kvec[r], Mlist[q][r][s], Mlist[qp][r][s]) -
                            expComp(llims[s], kvec[r], Mlist[q][r][s], Mlist[qp][r][s]);
                }//s
            }//r
            N += Cf[q]* std::conj(Cf[qp]) * (N1qq + N2qq + N3qq);
        }//qp
    }//q
    return N;
}


/**==============================__3-Vectors_of_Integers__===================================**/
/// construct the sets of allowed integers for a given fundamental domain with the option to
/// return each set of ns or the set of n-vectors in the "retNs" struct

eigenbasis_lin_T3::retNs eigenbasis_lin_T3::domainNvec(){
    // make lists of nx, ny, nz from -nMax to nMax
    std::vector< std::vector<int> > nvec(3);
    retNs Ns;
    for (int nx = -nMax; nx <= nMax; nx++){
        nvec[0].push_back(nx);
        for (int ny = -nMax; ny <= nMax; ny++){
            nvec[1].push_back(ny);
            for (int nz = -nMax; nz <= nMax; nz++){
                nvec[2].push_back(nz);
            }//nx
        }//ny
    }//nz
    Ns.ns = nvec;
    // do transpose, making list of the vectors (nx, ny, nz)
    std::vector< std::vector<int> > nvecTr;
    for (int nx : nvec[0]){
        std::vector<int> nv(3, 0);
        nv[0] = nx;
        for (int ny : nvec[1]){
            nv[1] = ny;
            for (int nz : nvec[2]){
                nv[2] = nz;
                nvecTr.push_back(nv);
            }//nz
        }//ny
    }//nx
    Ns.nsTr = nvecTr;
    return Ns;
}

/**==============================__Eigenmode_Computation__================================**/
/// computes the eigenmodes of a given fundamental domain and returns it and its components in the struct "ret"

/// inputs:
/// rT: the observer orientation as a unit vector
/// Rlss: the position vector magnitude, the radius of the last scattering surface

/// outputs:
/// ret: an instance of the struct retVals containing the k-vector, coefficients, and basis sum for each domain I

eigenbasis_lin_T3::retVals eigenbasis_lin_T3::eigBasis(double Rlss, std::vector<double> rT){
    /// terms needed to specify the eigenmode basis
    eigenbasis_lin_T3::retVals ret;                                     // struct of return values
    std::vector< std::complex<double> > Cf;                             // CEI output
    std::complex<double> nrmC;                                          // NConst output
    std::vector<double> kvec(3, 0);                                     // KEI output
    std::complex<double> esum = 0;                                      // eigenmode basis summation term
    std::vector< std::complex<double> > xiterm(6, 0);                   // plane wave terms of eigenmode summation
    std::vector< std::vector<double> > I(3, std::vector<double>(3, 0)); // Identity matrix for q=0 cases
    I[0][0] = 1; I[1][1] = 1; I[2][2] = 1;
    std::vector< std::vector< std::vector<double> > > Mlist;             //format M's as list to use as input
    std::vector< std::vector<int> > nvec = domainNvec().nsTr;
    int nsz = nvec.size();
    std::vector<double> tvec(3, 0);

//    E1
    if (domIndx == 1){
        Mlist.push_back(I);
        tvec[0] = 0; tvec[1] = 0; tvec[2] = lengths[2];
        Cf.push_back(1.0);
        for (int n = 0; n < nsz; n++){
            kvec = K(nvec[n], angles, lengths);
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, tvec, Mlist);
            esum = sum_EI(xiterm, Cf);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
        Mlist.clear();
    }
//    E2
    if (domIndx == 2){
        std::vector< std::complex<double> > Cf_non1;
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1.0; Mb[1][1] = -1.0; Mb[2][2] = 1.0;
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0] == 0 && nvec[n][1] == 0 && nvec[n][2]%2 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if ((nvec[n][0] > 0) || (nvec[n][0] == 0 && nvec[n][1] > 0)) {
                Mlist.push_back(I); Mlist.push_back(Mb);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E3
    if (domIndx == 3){
        std::vector< std::complex<double> > Cf_non1;
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][1] = 1; Mb[1][0] = -1; Mb[2][2] = 1;
        std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
        Mb2 = matPow(Mb, 2);
        std::vector< std::vector<double> > Mb3(3, std::vector<double>(3, 0));
        Mb3 = matPow(Mb, 3);
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0] == 0 && nvec[n][1] == 0 && nvec[n][2]%4 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] > 0 && nvec[n][1] >= 0) {
                Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2); Mlist.push_back(Mb3);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = 0.25*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E4
    if (domIndx == 4){
        std::vector< std::complex<double> > Cf_non1;
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = -0.5; Mb[2][2] = 1;
        std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
        Mb2 = matPow(Mb, 2);
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0] == 0 && nvec[n][1] == 0 && nvec[n][2]%3 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] > 0 && nvec[n][1] >= 0) {
                Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = lengths[2]/3.0;
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            }
            nrmC = NConst_hex(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E5
    if (domIndx == 5){
        std::vector< std::complex<double> > Cf_non1;
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = -0.5; Mb[2][2] = 1;
        std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
        Mb2 = matPow(Mb, 2);
        std::vector< std::vector<double> > Mb3(3, std::vector<double>(3, 0));
        Mb3 = matPow(Mb, 3);
        std::vector< std::vector<double> > Mb4(3, std::vector<double>(3, 0));
        Mb4 = matPow(Mb, 4);
        std::vector< std::vector<double> > Mb5(3, std::vector<double>(3, 0));
        Mb5 = matPow(Mb, 5);
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0] == 0 && nvec[n][1] == 0 && nvec[n][2]%6 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] > 0 && nvec[n][1] >= 0) {
                Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2); Mlist.push_back(Mb3); Mlist.push_back(Mb4); Mlist.push_back(Mb5);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = lengths[2]/6.0;
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            }
            nrmC = NConst_hex(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E6
    if (domIndx == 6){
        std::vector< std::complex<double> > Cf_non1;
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = -1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = -1;
        std::vector< std::vector<double> > Mc(3, std::vector<double>(3, 0));
        Mc[0][0] = -1; Mc[1][1] = -1; Mc[2][2] = 1;
        for (int n = 0; n < nsz; n++){
            if ( (nvec[n][0] > 0 && nvec[n][1] == 0 && nvec[n][2] == 0) ||
                 (nvec[n][1] > 0 && nvec[n][0] == 0 && nvec[n][2] == 0) ||
                 (nvec[n][2] > 0 && nvec[n][0] == 0 && nvec[n][1] == 0) ){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] > 0 && nvec[n][1] > 0) {
                Mlist.push_back(I); Mlist.push_back(Ma);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] =0 ;
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            } else if (nvec[n][0] == 0 && nvec[n][1] > 0 && nvec[n][2] > 0) {
                Mlist.push_back(I); Mlist.push_back(Mb);
                tvec[0] = 0; tvec[1] = 0.5*lengths[1]; tvec[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            } else if (nvec[n][0] > 0 && nvec[n][1] == 0 && nvec[n][2] > 0) {
                Mlist.push_back(I); Mlist.push_back(Mc);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0; tvec[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            }
            nrmC = NConst_hex(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E7
    if (domIndx == 7){
        std::vector< std::complex<double> > Cf_mata1, Cf_ta2;
        std::vector<double> tvec_a2(3, 0);
        std::vector< std::vector< std::vector<double> > > Mtemp;             //format M's as list to use as input
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0]%2 == 0 && nvec[n][1] == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][1] > 0) {
                Mlist.push_back(I); Mlist.push_back(Ma);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] = 0;
                tvec_a2[0] = 0.5*lengths[0]; tvec_a2[1] = -0.5*lengths[1]; tvec_a2[2] = 0;
                kvec = K(nvec[n], angles, lengths);
                Mtemp.push_back(Ma);
                Cf_mata1 = expMatArg(kvec, tvec, Mtemp);
                Mtemp.clear(); Mtemp.push_back(I);
                Cf_ta2 = expMatArg(kvec, tvec_a2, Mtemp);
                Mtemp.clear();
                Cf.push_back(Cf_mata1[1] * Cf_ta2[1]);
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E8
    if (domIndx == 8){
        std::vector< std::complex<double> > Cf_mata1, Cf_mata2, Cf_non1, Cf_mbtb;
        std::vector<double> tvec_a2(3, 0);
        std::vector<double> tvec_b(3, 0);
        std::vector< std::vector< std::vector<double> > > Mtemp;             //format M's as list to use as input
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = 1;
        std::vector< std::vector<double> > MaMb(3, std::vector<double>(3, 0));
        MaMb[0][0] = -1; MaMb[1][1] = -1; MaMb[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0]%2 == 0 && nvec[n][1] == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] == 0 && nvec[n][1] > 0 && nvec[n][2]%2 == 0) {
                Mlist.push_back(Ma);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] = 0;
                tvec_a2[0] = 0.5*lengths[0]; tvec_a2[1] = -0.5*lengths[0]; tvec_a2[2] = 0;
                kvec = K(nvec[n], angles, lengths);
                Cf_mata1 = expMatArg(kvec, tvec, Mlist);
                Cf_mata2 = expMatArg(kvec, tvec_a2, Mlist);
                Cf.push_back(Cf_mata1[1] * Cf_mata2[1]);
            } else if (nvec[n][0] > 0 && nvec[n][0]%2 == 0 && nvec[n][1] == 0) {
                Mlist.push_back(Mb);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            } else if (nvec[n][0] > 0 && nvec[n][1] > 0) {
                Mlist.push_back(Ma); Mlist.push_back(Mb); Mlist.push_back(MaMb);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] = 0;
                tvec_a2[0] = 0.5*lengths[0]; tvec_a2[1] = -0.5*lengths[0]; tvec_a2[2] = 0;
                tvec_b[0] = 0; tvec_b[1] = 0; tvec_b[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Mtemp.push_back(Ma);
                Cf_mata1 = expMatArg(kvec, tvec, Mtemp);
                Cf_mata2 = expMatArg(kvec, tvec_a2, Mtemp);
                Mtemp.clear(); Mtemp.push_back(Mb);
                Cf_mbtb = expMatArg(kvec, tvec_b, Mtemp);
                Mtemp.clear();
                Cf.push_back(Cf_mata1[1] * Cf_mata2[1] * Cf_mbtb[1]);
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E9
    if (domIndx == 9){
        std::vector< std::complex<double> > Cf_t1, Cf_t2, Cf_t3;
        std::vector<double> tvec_2(3, 0);
        std::vector<double> tvec_3(3, 0);
        std::vector< std::vector< std::vector<double> > > Mtemp;            //format M's as list to use as input
        Mtemp.push_back(I);
        Cf.push_back(1.0);
        std::vector< std::vector<double> > M(3, std::vector<double>(3, 0));
        M[0][0] = 1; M[1][1] = -1; M[2][2] = 1;
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0]%2 == 0 && nvec[n][1] == 0 && nvec[n][2]%2 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][1] > 0 && (nvec[n][0] + nvec[n][1])%2 == nvec[n][2]%2) {
                Mlist.push_back(M);
                kvec = K(nvec[n], angles, lengths);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[0]; tvec[2] = 0;
                tvec_2[0] = 0.5*lengths[0]; tvec_2[1] = -0.5*lengths[0]; tvec_2[2] = 0;
                tvec_3[0] = 0; tvec_3[1] = 0; tvec_3[2] = lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_t1 = expMatArg(kvec, tvec, Mlist);
                Cf_t2 = expMatArg(kvec, tvec_2, Mlist);
                Cf_t3 = expMatArg(kvec, tvec_3, Mlist);
                Cf.push_back(Cf_t1[1] * Cf_t1[2] * Cf_t3[1]);
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
//    E10
    if (domIndx == 10){
        std::vector< std::complex<double> > Cf_mata1, Cf_mata2, Cf_non1, Cf_mbtb;
        std::vector<double> tvec_a2(3, 0);
        std::vector<double> tvec_b(3, 0);
        std::vector< std::vector< std::vector<double> > > Mtemp;             //format M's as list to use as input
        Cf.push_back(1.0);
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = -1; Mb[2][2] = 1;
        std::vector< std::vector<double> > MaMb(3, std::vector<double>(3, 0)); // M_A . M_B
        MaMb[0][0] = -1; MaMb[1][1] = 1; MaMb[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        for (int n = 0; n < nsz; n++){
            if (nvec[n][0] == 0 && nvec[n][1] == 0 && nvec[n][2]%2 == 0){
                Mlist.push_back(I);
                kvec = K(nvec[n], angles, lengths);
            } else if (nvec[n][0] == 0 && nvec[n][1] > 0 && nvec[n][2]%2 == 0 && nvec[n][1]%2 == nvec[n][2]%2) {
                Mlist.push_back(Ma);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] = 0;
                tvec_a2[0] = 0.5*lengths[0]; tvec_a2[1] = -0.5*lengths[0]; tvec_a2[2] = 0;
                kvec = K(nvec[n], angles, lengths);
                Cf_mata1 = expMatArg(kvec, tvec, Mlist);
                Cf_mata2 = expMatArg(kvec, tvec_a2, Mlist);
                Cf.push_back(Cf_mata1[1] * Cf_mata2[1]);
            } else if (nvec[n][0] > 0 && nvec[n][0]%2 == 0 && nvec[n][1] == 0) {
                Mlist.push_back(Mb);
                tvec[0] = 0; tvec[1] = 0; tvec[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Cf_non1 = expMatArg(kvec, tvec, Mlist);
                for (int ii; ii < Cf_non1.size(); ii++){
                    Cf.push_back(Cf_non1[ii]);
                }
            } else if (nvec[n][0] == 0 && nvec[n][1] > 0) {
                Mlist.push_back(Ma); Mlist.push_back(Mb); Mlist.push_back(MaMb);
                tvec[0] = 0.5*lengths[0]; tvec[1] = 0.5*lengths[1]; tvec[2] = 0;
                tvec_a2[0] = 0.5*lengths[0]; tvec_a2[1] = -0.5*lengths[0]; tvec_a2[2] = 0;
                tvec_b[0] = 0; tvec_b[1] = 0; tvec_b[2] = 0.5*lengths[2];
                kvec = K(nvec[n], angles, lengths);
                Mtemp.push_back(Ma);
                Cf_mata1 = expMatArg(kvec, tvec, Mtemp);
                Cf_mata2 = expMatArg(kvec, tvec_a2, Mtemp);
                Mtemp.clear(); Mtemp.push_back(Mb);
                Cf_mbtb = expMatArg(kvec, tvec_b, Mtemp);
                Mtemp.clear();
                Cf.push_back(Cf_mata1[1] * Cf_mata2[1] * Cf_mbtb[1]);
            }
            nrmC = NConst_rprism(lengths, kvec, Cf, Mlist);
            xiterm = expMatArg(kvec, rT, Mlist);
            Mlist.clear();
            esum = sum_EI(Cf, xiterm);
            ret.k.push_back(kvec);
            ret.basis.push_back(nrmC*std::conj(nrmC)*esum);
        }//n
    }
    return ret;
}

/// below are functions for setting and getting the private variables of the eigenbasis_lin_T3 class

void eigenbasis_lin_T3::set_nMax(int nMaxIn){
    nMax = nMaxIn;
}

int eigenbasis_lin_T3::get_nMax(){
    return nMax;
}

void eigenbasis_lin_T3::set_domIndx(int domIndxIn){
    domIndx = domIndxIn;
}

int eigenbasis_lin_T3::get_domIndx(){
    return domIndx;
}

void eigenbasis_lin_T3::set_angles(std::vector<double> anglesIn){
    angles[0] = anglesIn[0];
    angles[1] = anglesIn[1];
    angles[2] = anglesIn[2];
}

std::vector<double> eigenbasis_lin_T3::get_angles(){
    return angles;
}

void eigenbasis_lin_T3::set_lengths(std::vector<double> lengthsIn){
    lengths[0] = lengthsIn[0];
    lengths[1] = lengthsIn[1];
    lengths[2] = lengthsIn[2];
}

std::vector<double> eigenbasis_lin_T3::get_lengths(){
    return lengths;
}
