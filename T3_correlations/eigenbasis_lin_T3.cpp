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

/// parallelepiped k
std::vector<double> Kll(std::vector<int> nvec, std::vector<double> angles, std::vector<double> lengths){
    std::vector<double> k(3, 0);
    k[0] = 2 * M_PI * nvec[0] / lengths[0];
    k[1] = 2 * M_PI * ((nvec[1] / (std::sin(angles[0]) * lengths[1])) - (nvec[0] / (std::tan(angles[0]) * lengths[0])));
    k[2] = 2 * M_PI * ((lengths[2] * (lengths[0] * nvec[1] * std::cos(angles[1]) + lengths[1] * nvec[0] * std::cos(angles[2])) / std::tan(angles[0])) -
                       (lengths[2] * (lengths[1] * nvec[0] * std::cos(angles[1]) + lengths[0] * nvec[1] * std::cos(angles[2])) / std::sin(angles[0])) +
                       lengths[0] * lengths[1] * nvec[2] * std::sin(angles[0]) ) /
           (lengths[0] * lengths[1] * lengths[2] * std::sqrt(1 + 2 * std::cos(angles[0])*std::cos(angles[1])*std::cos(angles[2]) - std::cos(angles[0])*std::cos(angles[0]) - std::cos(angles[1])*std::cos(angles[1]) - std::cos(angles[2])*std::cos(angles[2])) );
    return k;
}

/// regular hexagonal prism k
std::vector<double> Khx(std::vector<int> nvec, std::vector<double> angles, std::vector<double> lengths){
    std::vector<double> k(3, 0);
    double L1 = -1 / (std::sqrt(1 + (lengths[0]/lengths[1])*(lengths[0]/lengths[1])));
    double L2 = 1 / (lengths[1]*std::sqrt(1 + (lengths[0]/lengths[1])*(lengths[0]/lengths[1])));
    k[0] = 2 * M_PI * nvec[0] / lengths[0];
    k[1] = 2 * M_PI * nvec[1] / L1;
    k[2] = 2 * M_PI * (L2*lengths[2]*nvec[0]*std::cos(angles[1]) + L1*lengths[2]*nvec[1]*std::cos(angles[2]) - L1*L2*nvec[2]) / (L1*L2*lengths[2] * std::sqrt(std::sin(angles[1])*std::sin(angles[1]) - std::cos(angles[2])*std::cos(angles[2]) ));
    return k;
}

/// rhombic dodecahedron k
//std::vector<double> Krd(std::vector<int> nvec, std::vector<double> angles, std::vector<double> lengths){
//    std::vector<double> k(3, 0);
//    k[0] = 2 * M_PI * nvec[0] / lengths[0];
//    k[1] = 2 * M_PI * nvec[1] / lengths[1];
//    k[2] = 2 * M_PI * nvec[2] / lengths[2];
//    return k;
//}

/**================================__Coefficients__=================================**/
/// some notation:
/// CEI -> Ith fundamental domain's eigenmode's basis coefficients

/// inputs:
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI
/// nvec: 3-vector of elements of each allowed set of integers

/// outputs:
/// coeffVec: complex-valued 3-vector containing each a_q needed for representing
///           the eigenmodes as linear combinations of the covering space eigenmodes of domain I

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

std::vector< std::complex<double> > CE1(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;

    std::complex<double> cterm = 1.0;
    coeffVec.push_back(cterm);

    return coeffVec;
}

std::vector< std::complex<double> > CE2(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1.0; Mb[1][1] = -1.0; Mb[2][2] = 1.0;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm = 0;
    std::complex<double> epow = 0;
    std::vector<double> t(3, 0);
    t[0] = 0; t[1] = 0; t[2] = 0.5*lengths[2];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
             epow += kvec[i]*Mb[i][j]*t[j]; //a_1
        }
    }
    cterm = std::exp(imu*epow);
    coeffVec.push_back(cterm);

    return coeffVec;
}

std::vector< std::complex<double> > CE3(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][1] = 1; Mb[1][1] = -1; Mb[2][2] = 1;
    std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
    Mb2 = matPow(Mb, 2);
    std::vector< std::vector<double> > Mb3(3, std::vector<double>(3, 0));
    Mb3 = matPow(Mb, 3);
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3;
    std::complex<double> epow1, epow2, epow3;
    std::vector<double> t(3, 0);
    t[0] = 0; t[1] = 0; t[2] = 0.25*lengths[2];
    cterm1 = 0;cterm2 = 0; cterm3 = 0;
    epow1 = 0; epow2 = 0; epow3 = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1 += kvec[i]*Mb[i][j]*t[j]; //a_1
            epow2 += kvec[i]*Mb2[i][j]*t[j]; //a_2
            epow3 += kvec[i]*Mb3[i][j]*t[j]; //a_3
        }
    }
    cterm1 = std::exp(imu*epow1); //a_1
    cterm2 = std::exp(imu*epow2); //a_2
    cterm3 = std::exp(imu*epow3); //a_3
    cterm2 *= cterm1;
    cterm3 *= cterm2;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);
    coeffVec.push_back(cterm3);

    return coeffVec;
}

std::vector< std::complex<double> > CE4(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = -0.5; Mb[2][2] = 1;
    std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
    Mb2 = matPow(Mb, 2);
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3;
    std::complex<double> epow1, epow2, epow3;
    std::vector<double> t(3, 0);
    t[0] = 0; t[1] = 0; t[2] = lengths[2]/3;
    epow1 = 0; epow2 = 0; epow3 = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1 += kvec[i]*t[j];
            epow2 += kvec[i]*Mb[i][j]*t[j];//a_1
            epow3 += kvec[i]*Mb2[i][j]*t[j]; //a_2
        }
    }
    cterm1 = std::exp(imu*epow1);
    cterm2 = std::exp(imu*epow2);
    cterm3 = std::exp(imu*epow3);
    cterm2 *= cterm1;
    cterm3 *= cterm2;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);

    return coeffVec;
}

std::vector< std::complex<double> > CE5(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = 0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = 0.5; Mb[2][2] = 1;
    std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
    Mb2 = matPow(Mb, 2);
    std::vector< std::vector<double> > Mb3(3, std::vector<double>(3, 0));
    Mb3 = matPow(Mb, 3);
    std::vector< std::vector<double> > Mb4(3, std::vector<double>(3, 0));
    Mb4 = matPow(Mb, 4);
    std::vector< std::vector<double> > Mb5(3, std::vector<double>(3, 0));
    Mb5 = matPow(Mb, 5);
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3, cterm4, cterm5, cterm6;
    std::complex<double> epow1, epow2, epow3, epow4, epow5, epow6;
    std::vector<double> t(3, 0);
    t[0] = 0; t[1] = 0; t[2] = lengths[2]/6;
    epow1 = 0; epow2 = 0; epow3 = 0; epow4 = 0; epow5 = 0; epow6 = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1 += kvec[i]*t[j];
            epow2 += kvec[i]*Mb[i][j]*t[j]; //a_1
            epow3 += kvec[i]*Mb2[i][j]*t[j]; //a_2
            epow4 += kvec[i]*Mb3[i][j]*t[j]; //a_3
            epow5 += kvec[i]*Mb4[i][j]*t[j]; //a_4
            epow6 += kvec[i]*Mb5[i][j]*t[j]; //a_5
        }
    }
    cterm1 = std::exp(imu*epow1);
    cterm2 = std::exp(imu*epow2);
    cterm3 = std::exp(imu*epow3);
    cterm4 = std::exp(imu*epow4);
    cterm5 = std::exp(imu*epow5);
    cterm6 = std::exp(imu*epow6);
    cterm2 *= cterm1;
    cterm3 *= cterm2;
    cterm4 *= cterm3;
    cterm5 *= cterm4;
    cterm6 *= cterm5;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);
    coeffVec.push_back(cterm3);
    coeffVec.push_back(cterm4);
    coeffVec.push_back(cterm5);
    coeffVec.push_back(cterm6);

    return coeffVec;
}

std::vector< std::complex<double> > CE6(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = -1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = -1;
    std::vector< std::vector<double> > Mc(3, std::vector<double>(3, 0));
    Mc[0][0] = -1; Mc[1][1] = -1; Mc[2][2] = 1;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3, cterm4;
    std::complex<double> epow1, epow2, epow3, epow4;
    std::vector<double> ta(3, 0), tb(3, 0), tc(3, 0);
    ta[0] = 0.5*lengths[0]; ta[1] = 0.5*lengths[1]; ta[2] = 0;              //ta
    tb[0] = 0;              tb[1] = 0.5*lengths[1]; tb[2] = 0.5*lengths[2]; //tb
    tc[0] = 0.5*lengths[0]; tc[1] = 0;              tc[2] = 0.5*lengths[2]; //tc
    epow1 = 0; epow2 = 0; epow3 = 0; epow4 = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1 += kvec[i]*ta[j]; //a_1
            epow2 += kvec[i]*Ma[i][j]*ta[j]; //a_1
            epow3 += kvec[i]*Mb[i][j]*tb[j]; //a_2
            epow4 += kvec[i]*Mc[i][j]*tc[j]; //a_3
        }
    }
    cterm1 = std::exp(imu*epow1);
    cterm2 = std::exp(imu*(epow2 - epow3));
    cterm3 = std::exp(imu*(epow3 - epow4));
    cterm4 = std::exp(imu*(epow4 - epow2));
    cterm2 *= cterm1;
    cterm3 *= cterm2;
    cterm4 *= cterm3;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);
    coeffVec.push_back(cterm3);
    coeffVec.push_back(cterm4);
    return coeffVec;
}

std::vector< std::complex<double> > CE7(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2;
    std::complex<double> epow1x, epow2x;
    std::complex<double> epow1y, epow2y;
    std::vector<double> ta1(3, 0), ta2(3, 0);
    ta1[0] = 0.5*lengths[0]; ta1[1] = 0.5*lengths[1]; ta1[2] = 0;
    ta2[0] = 0.5*lengths[0]; ta2[1] = -0.5*lengths[1]; ta2[2] = 0;
    epow1x = 0; epow2x = 0; epow1y = 0; epow2y = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1x += kvec[i]*ta1[j];
            epow1y += kvec[i]*ta2[j];
            epow2x += kvec[i]*Ma[i][j]*ta1[j]; //a_1
            epow2y += kvec[i]*Ma[i][j]*ta2[j]; //a_1
        }
    }
    cterm1 = std::exp(imu*(epow1x + epow1y));
    cterm2 = std::exp(imu*(epow2x + epow2y));
    cterm2 *= cterm1;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);

    return coeffVec;
}

std::vector< std::complex<double> > CE8(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = 1;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3, cterm4;
    std::complex<double> epow1x, epow2x;
    std::complex<double> epow1y, epow2y;
    std::complex<double> epow1z, epow2z;
    std::vector<double> ta1(3, 0), ta2(3, 0), tb(3, 0);
    ta1[0] = 0.5*lengths[0]; ta1[1] = 0.5*lengths[1]; ta1[2] = 0;
    ta2[0] = 0.5*lengths[0]; ta2[1] = -0.5*lengths[1]; ta2[2] = 0;
    tb[0] = 0.5*lengths[0]; tb[1] = -0.5*lengths[1]; tb[2] = 0;
    epow1x = 0; epow2x = 0; epow1y = 0; epow2y = 0; epow1z = 0; epow2z = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1x += kvec[i]*ta1[j];
            epow1y += kvec[i]*ta2[j];
            epow1z += kvec[i]*tb[j];
            epow2x += kvec[i]*Ma[i][j]*ta1[j]; //a_1
            epow2y += kvec[i]*Ma[i][j]*ta2[j]; //a_1
            epow2z += kvec[i]*Mb[i][j]*tb[j];
        }
    }
    cterm1 = std::exp(imu*(epow1x + epow1y + epow1z));
    cterm2 = std::exp(imu*(epow2x + epow2y));
    cterm3 = std::exp(imu*epow2z);
    cterm4 = std::exp(imu*(epow2x + epow2y + epow2z));
//    cterm2 *= cterm1;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);
    coeffVec.push_back(cterm3);
    coeffVec.push_back(cterm4);

    return coeffVec;
}

std::vector< std::complex<double> > CE9(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > M(3, std::vector<double>(3, 0));
    M[0][0] = 1; M[1][1] = -1; M[2][2] = 1;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2;
    std::complex<double> epow1x, epow2x;
    std::complex<double> epow1y, epow2y;
    std::vector<double> ta1(3, 0), ta2(3, 0);
    ta1[0] = 0.5*lengths[0]; ta1[1] = 0.5*lengths[1]; ta1[2] = 0;
    ta2[0] = 0.5*lengths[0]; ta2[1] = -0.5*lengths[1]; ta2[2] = 0;
    epow1x = 0; epow2x = 0; epow1y = 0; epow2y = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1x += kvec[i]*ta1[j];
            epow1y += kvec[i]*ta2[j];
            epow2x += kvec[i]*M[i][j]*ta1[j]; //a_1
            epow2y += kvec[i]*M[i][j]*ta2[j]; //a_1
        }
    }
    cterm1 = std::exp(imu*(epow1x + epow1y));
    cterm2 = std::exp(imu*(epow2x + epow2y));
    cterm2 *= cterm1;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);

    return coeffVec;
}

std::vector< std::complex<double> > CE10(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths){
    std::vector< std::complex<double> >  coeffVec;
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = -1; Mb[2][2] = 1;
//    coeffVec.push_back(1.0); //a_0
    std::complex<double> cterm1, cterm2, cterm3, cterm4;
    std::complex<double> epow1x, epow2x;
    std::complex<double> epow1y, epow2y;
    std::complex<double> epow1z, epow2z;
    std::vector<double> ta1(3, 0), ta2(3, 0), tb(3, 0);
    ta1[0] = 0.5*lengths[0]; ta1[1] = 0.5*lengths[1]; ta1[2] = 0;
    ta2[0] = 0.5*lengths[0]; ta2[1] = -0.5*lengths[1]; ta2[2] = 0;
    tb[0] = 0.5*lengths[0]; tb[1] = -0.5*lengths[1]; tb[2] = 0;
    epow1x = 0; epow2x = 0; epow1y = 0; epow2y = 0; epow1z = 0; epow2z = 0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            epow1x += kvec[i]*ta1[j];
            epow1y += kvec[i]*ta2[j];
            epow1z += kvec[i]*tb[j];
            epow2x += kvec[i]*Ma[i][j]*ta1[j]; //a_1
            epow2y += kvec[i]*Ma[i][j]*ta2[j]; //a_1
            epow2z += kvec[i]*Mb[i][j]*tb[j];
        }
    }
    cterm1 = std::exp(imu*(epow1x + epow1y + epow1z));
    cterm2 = std::exp(imu*(epow2x + epow2y));
    cterm3 = std::exp(imu*epow2z);
    cterm4 = std::exp(imu*(epow2x + epow2y + epow2z));
//    cterm2 *= cterm1;
    coeffVec.push_back(cterm1);
    coeffVec.push_back(cterm2);
    coeffVec.push_back(cterm3);
    coeffVec.push_back(cterm4);

    return coeffVec;
}

/**==============================__Normalization_Constants__===================================**/
/// inputs:
/// lengths: 3 vector of fundamental domain side lengths
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI; must match domIndx to NEI as well
/// coeffs: coefficients of the eigenmode's basis; outputs of CEI functions
/// transMats: vector of transformation matrices or products of transfomation matrices

/// output:
/// 1.0 / sqrt(N): The normalization constant for the linear representation of eigenmode functions for the domain I

std::complex<double> NConst (std::vector<double> lengths, std::vector<double> kvec, std::vector< std::complex<double> > coeffs, std::vector< std::vector< std::vector<double> > > transMats){
    int q_max = transMats.size();
    std::complex<double> N = 0.0; // technically 1/normconst^2; evaluate at the end
    std::complex<double> imPart = 0.0; //imaginary part of normalization constant
    std::vector< std::vector<double> > v(q_max, std::vector<double>(3, 0)); // k.M_q; vector of vectors so we can select M_q
    std::vector< std::complex<double> > expTrm(q_max, 1.0); //e^iv_q.L
    std::complex<double> expqp = 1.0; //e^iv_q'.L
    double denom = 1.0; // 1/prod_p(v_qp - v_q'p)
    N = q_max * lengths[0] * lengths[1] * lengths[2];
    /// compute k.M_q and k.M_q'
    for (int a = 0; a < 3; ++a){
        for (int b = 0; b < 3; ++b){
            for (int q = 0; q < q_max; ++q){
                v[q][a] += kvec[a] * transMats[q][a][b];
            }//q
        }//b
    }//a
    /// compute 1 / /prod_p(v_{qp} - v_{q'p})
    for (int a = 0; a < 3; ++a){
        for (int q = 0; q < q_max; ++q){
            for (int qp = 0; qp < q_max; ++qp){
                denom *= 1 / (v[q][a] - v[qp][a]);
            }//qp
        }//q
    }//a
    /// compute e^iv_q.L and e^iv_q'.L
    for (int q = 0; q < q_max; ++q){
        for (int a = 0; a < 3; ++a){
            expTrm[q] *= std::exp(-imu*v[q][a]*lengths[a]);
        }//a
    }//q
    /// compute /prod_p(e^iv_{qp}L_p - e^iv_{q'p}L_p}) and /prod_p(e^iv_{q'p}L_p - e^iv_{qp}L_p})
    std::complex<double> prod1 = 1.0;
    std::complex<double> prod2 = 1.0;
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            if (q != qp){
                for (int p = 0; p < 3; ++p){
                    prod1 *= (std::exp(imu*v[q][p]*lengths[p]) - std::exp(imu*v[qp][p]*lengths[p]));
                    prod2 *= (std::exp(imu*v[qp][p]*lengths[p]) - std::exp(imu*v[q][p]*lengths[p]));
                }//p
                imPart += coeffs[q]*conj(coeffs[qp])*(expTrm[q]*prod1 + expTrm[qp]*prod2);
            }//if
            else imPart = 0;
        }//qp
    }//q
    imPart *= denom;
    N += imPart;
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}
/**==============================__3-Vectors_of_Integers__===================================**/
/// construct the sets of allowed integers for a given fundamental domain with the option to
/// return each set of ns or the set of n-vectors in the "retNs" struct

eigenbasis_lin_T3::retNs eigenbasis_lin_T3::domainNvec(){
    std::vector< std::vector<int> > nvec(3);
    retNs Ns;
    if (domIndx == 1){
        for (int n = -nMax; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
            nvec[2].push_back(n);
        }
    }
    else if (domIndx == 2){
        nvec[0].push_back(0);
        for (int n = 1; n <= nMax; n++) nvec[1].push_back(n);
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 3){
        nvec[1].push_back(0);
        for (int n = 1; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
        }
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 4){
        nvec[1].push_back(0);
        for (int n = 1; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
        }
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 5){
        nvec[1].push_back(0);
        for (int n = 1; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
        }
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 6){
        for (int n = 1; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
        }
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 7){
        for (int n = 1; n <= nMax; n++) nvec[1].push_back(n);
        for (int n = -nMax; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[2].push_back(n);
        }
    }
    else if (domIndx == 8){
        for (int n = 1; n <= nMax; n++){
            nvec[0].push_back(n);
            nvec[1].push_back(n);
        }
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else if (domIndx == 9){
        for (int ny = 1; ny <= nMax; ny++) nvec[1].push_back(ny);
        for (int nx = -nMax; nx <= nMax; nx++) nvec[0].push_back(nx);
        for (int ny = 1; ny <= nMax; ny++){
            for (int nx = -nMax; nx <= nMax; nx++){
                for (int nz = -nMax; nz <= nMax; nz++){
                    if ( (nx+ny)%2 == nz%2 ){ // (nx+ny)mod2 def nzmod2
                        nvec[2].push_back(nz);
                    }
                }
            }
        }
    }
    else if (domIndx == 10){
        nvec[0].push_back(0);
        for (int n = 1; n <= nMax; n++) nvec[1].push_back(n);
        for (int n = -nMax; n <= nMax; n++) nvec[2].push_back(n);
    }
    else {std::cout << "error: domainNvec() only handles E1-E10 quotient spaces";}
    Ns.ns = nvec;
    // do transpose
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
    eigenbasis_lin_T3::retVals ret;         // struct of return values
    std::vector< std::complex<double> > Cf; // CEI output
    std::complex<double> nrmC;              // NConst output
    std::vector<double> kvec(3, 0);         // KEI output
    std::vector< std::vector<int> > nvec;   // array of domainNvec().nsTr outputs
    std::complex<double> esum = 0;          // eigenmode basis summation term
    std::vector< std::complex<double> > xiterm(6, 0); // plane wave terms of eigenmode summation
    std::vector< std::vector<double> > I(3, std::vector<double>(3, 0)); // Identity matrix for q=0 cases
    I[0][0] = 1; I[1][1] = 1; I[2][2] = 1;

//    E1
    if (domIndx == 1){
        esum = 0;
        nvec = domainNvec().nsTr;
        std::vector< std::vector<double> > M (3, std::vector<double>(3, 0));
        M[0][0] = 1; M[1][1] = 1; M[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(M);
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE1(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*M[j][k]*rT[k]);
                    }//k
                }//j
                esum += Cf[i]*xiterm[0];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E2
    if (domIndx == 2){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1.0; Mb[1][1] = -1.0; Mb[2][2] = 1.0;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Mb);
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE2(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E3
    if (domIndx == 3){
        esum = 0;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][1] = 1; Mb[1][0] = -1; Mb[2][2] = 1;
        std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
        Mb2 = matPow(Mb, 2);
        std::vector< std::vector<double> > Mb3(3, std::vector<double>(3, 0));
        Mb3 = matPow(Mb, 3);
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2); Mlist.push_back(Mb3);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE3(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                xiterm[3] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb2[j][k]*rT[k]);
                        xiterm[3] += std::exp(imu*kvec[j]*Mb3[j][k]*rT[k]);

                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E4
    if (domIndx == 4){
        esum = 0;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = -0.5; Mb[2][2] = 1;
        std::vector< std::vector<double> > Mb2(3, std::vector<double>(3, 0));
        Mb2 = matPow(Mb, 2);
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Khx(nvec[n], angles, lengths);
            Cf = CE4(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb2[j][k]*rT[k]);
                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E5
    if (domIndx == 5){
        esum = 0;
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
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Mb); Mlist.push_back(Mb2); Mlist.push_back(Mb3); Mlist.push_back(Mb4); Mlist.push_back(Mb5);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Khx(nvec[n], angles, lengths);
            Cf = CE5(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                xiterm[3] = 0;
                xiterm[4] = 0;
                xiterm[5] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb2[j][k]*rT[k]);
                        xiterm[3] += std::exp(imu*kvec[j]*Mb3[j][k]*rT[k]);
                        xiterm[4] += std::exp(imu*kvec[j]*Mb4[j][k]*rT[k]);
                        xiterm[5] += std::exp(imu*kvec[j]*Mb5[j][k]*rT[k]);

                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E6
    if (domIndx == 6){
        esum = 0;
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = -1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = -1;
        std::vector< std::vector<double> > Mc(3, std::vector<double>(3, 0));
        Mc[0][0] = -1; Mc[1][1] = -1; Mc[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Ma); Mlist.push_back(Mb); Mlist.push_back(Mc);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE6(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                xiterm[3] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Ma[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[3] += std::exp(imu*kvec[j]*Mc[j][k]*rT[k]);

                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E7
    if (domIndx == 7){
        esum = 0;
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        nvec = domainNvec().nsTr;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Ma);
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE7(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Ma[j][k]*rT[k]);
                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(esum);
    }
//    E8
    if (domIndx == 8){
        esum = 0;
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = 1;
        std::vector< std::vector<double> > MaMb(3, std::vector<double>(3, 0));
        MaMb[0][0] = -1; MaMb[1][1] = -1; MaMb[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Ma); Mlist.push_back(Mb); Mlist.push_back(MaMb);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE8(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                xiterm[3] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Ma[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[3] += std::exp(imu*kvec[j]*MaMb[j][k]*rT[k]);

                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(esum);
    }
//    E9
    if (domIndx == 9){
        esum = 0;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = 1; Mb[1][1] = -1; Mb[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Mb);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE9(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(esum);
    }
//    E10
    if (domIndx == 10){
        esum = 0;
        std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
        Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
        std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
        Mb[0][0] = -1; Mb[1][1] = -1; Mb[2][2] = 1;
        std::vector< std::vector<double> > MaMb(3, std::vector<double>(3, 0)); // M_A . M_B
        MaMb[0][0] = -1; MaMb[1][1] = 1; MaMb[2][2] = 1;
        std::vector< std::vector< std::vector<double> > > Mlist; //format M as list for input to NConst
        Mlist.push_back(I); Mlist.push_back(Ma); Mlist.push_back(Mb); Mlist.push_back(MaMb);
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE10(kvec, nvec[n], lengths);
            nrmC = NConst(lengths, kvec, Cf, Mlist);
            for (int i = 0; i < Cf.size(); i++){
                xiterm[0] = 0;
                xiterm[1] = 0;
                xiterm[2] = 0;
                xiterm[3] = 0;
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        xiterm[0] += std::exp(imu*kvec[j]*rT[k]);
                        xiterm[1] += std::exp(imu*kvec[j]*Ma[j][k]*rT[k]);
                        xiterm[2] += std::exp(imu*kvec[j]*Mb[j][k]*rT[k]);
                        xiterm[3] += std::exp(imu*kvec[j]*MaMb[j][k]*rT[k]);

                    }//k
                }//j
                esum += Cf[i]*xiterm[i];
            }//i
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
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

