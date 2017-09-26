#include "EigBase.hpp"
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
///// condition: index of the periodicity condition on the set of integers allowed for each wave-vector

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
    k[0] = 2 * M_PI * nvec[0] / lengths[0];
    k[1] = 2 * M_PI * (nvec[0] + 2 * nvec[1]) / (std::sqrt(3) * lengths[0]);
    k[2] = 2 * M_PI * (lengths[0 * nvec[2] - lengths[2] * nvec[0] * std::cos(angles[1])) / (lengths[0] * lengths[2] * std::sin(angles[1]));
    return k;
}

/// rhombic dodecahedron k
std::vector<double> Krd(std::vector<int> nvec, std::vector<double> angles, std::vector<double> lengths){
    std::vector<double> k(3, 0);
    k[0] = 2 * M_PI * nvec[0] / lengths[0];
    k[1] = 2 * M_PI * nvec[1] / lengths[1];
    k[2] = 2 * M_PI * nvec[2] / lengths[2];
    return k;
}

/**================================__Coefficients__=================================**/
/// some notation:
/// CEI -> Ith fundamental domain's eigenmode's basis coefficients

/// inputs:
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI
/// nvec: 3-vector of elements of each allowed set of integers
/// lMax: maximum azimuthal number of the spherical harmonics

double norm(std::vector<double> vec){ //compute 3-vector L2-norm
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

std::vector<double> hat(std::vector<double> vec){ //compute unit 3-vector
    std::vector<double> vhat(3, 0);
    if (vec[0] == 0 && vec[1] == 0 && vec[2] == 0){ //avoid x/0; if no vector then no direction so output zero vector
        vhat[0] = 0; vhat[1] = 0; vhat[2] = 0;
    }
    else {
        for (int i = 0; i < 3; i++){
            vhat[i] = (1 / norm(vec)) * vec[i];
        }
    }
    return vhat;
}

double polarAngle(std::vector<double> coordVec){ //compute polar angle of 3-vector in Cartesian coordinates
    if (coordVec[0] == 0 && coordVec[1] == 0 && coordVec[2] == 0){ //avoid x/0
        return 0;
    }
    else return acos( coordVec[2] / sqrt(coordVec[0]*coordVec[0] + coordVec[1]*coordVec[1] + coordVec[2]*coordVec[2]) );
}

double azimuthalAngle(std::vector<double> coordVec){ //compute azimuthal angle of a 3-vector in Cartesian coordinates
    if (coordVec[0] == 0 || coordVec[1] == 0){ //avoid x/0
        return 0;
    }
    else return tan( coordVec[1] / coordVec[0] );
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

std::vector< std::complex<double> > CE1(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    std::complex<double> YlmStar = 0; //spherical harmonic complex conjugate
//    int N = 1;
    double polangle = polarAngle(khat);
    double azimangle = azimuthalAngle(khat);
    for (int l = 2; l <= lMax; ++l){
        for (int m = -l; m <= l; ++m){
            YlmStar = std::conj(Y_lm(l, m, polangle, azimangle));
            coeffVec.push_back(std::pow(imu, l) * YlmStar);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE2(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    std::vector< std::complex<double> > YlmStar(2,0);
    std::complex<double> coeffTerm = 0;
    int N = 2;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1.0; Mb[1][1] = -1.0; Mb[2][2] = 1.0;
    std::vector<double> k_Mb(3, 0);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            k_Mb[i] += khat[j] * Mb[j][i];
        }
    }
    polangle[1] = polarAngle(k_Mb);
    azimangle[1]= azimuthalAngle(k_Mb);
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            YlmStar[0] = std::conj(Y_lm(l, m, polangle[0], azimangle[0]));
            YlmStar[1] = std::conj(Y_lm(l, m, polangle[1], azimangle[1]));
            coeffTerm = std::pow(imu, l) * ( YlmStar[0] - YlmStar[1]);
            coeffVec.push_back(coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE3(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 4;
    std::vector< std::complex<double> > YlmStar(N,0);
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Mb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][1] = 1; Mb[1][1] = -1; Mb[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Mb[n]);
        azimangle[n] = azimuthalAngle(k_Mb[n]);
    }
    for (int l = 2; l <= lMax; ++l){
        for (int m = -l; m < l + 1; ++m){
            for (int j = 0; j < N; ++j){
                YlmStar[j] = std::conj(Y_lm(l, m, polangle[j], azimangle[j]));
                coeffTerm += std::pow(imu,j*nvec[2]) * YlmStar[j];
            }
            coeffVec.push_back(std::pow(imu, l) * coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE4(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 2;
    std::vector< std::complex<double> > YlmStar(N);
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Mb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = -0.5; Mb[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Mb[n]);
        azimangle[n] = azimuthalAngle(k_Mb[n]);
    }
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            for (int j = 0; j < N; ++j){
                YlmStar[j] = conj(Y_lm(l, m, polangle[j], azimangle[j]));
                coeffTerm += std::pow(std::exp(imu * (2*M_PI/3) ), j*nvec[2]) * YlmStar[j];
            }
            coeffVec.push_back(std::pow(imu, l) * coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE5(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 2;
    std::vector< std::complex<double> > YlmStar(N);
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Mb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = 0.5; Mb[0][1] = -sqrt(3)/2; Mb[1][0] = sqrt(3)/2; Mb[1][1] = 0.5; Mb[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Mb[n]);
        azimangle[n] = azimuthalAngle(k_Mb[n]);
    }
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            for (int j = 0; j < N; ++j){
                YlmStar[j] = std::conj(Y_lm(l, m, polangle[j], azimangle[j]));
                coeffTerm += std::pow(std::exp(imu * (2*M_PI / 6)),j*nvec[2]) * YlmStar[j];
            }
            coeffVec.push_back(std::pow(imu, l) * coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE6(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 2;
    std::vector< std::complex<double> > YlmStar(4); //4 instead of N here.we use the same container differently later
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(4), azimangle(4);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Ma(N, std::vector<double>(3, 0)), k_Mb(N, std::vector<double>(3, 0)), k_Mc(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = -1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = -1;
    std::vector< std::vector<double> > Mc(3, std::vector<double>(3, 0));
    Mc[0][0] = -1; Mc[1][1] = -1; Mc[2][2] = 1;
    double polA, polB, polC, azimA, azimB, azimC;
    std::vector<double> kma(3, 0), kmb(3, 0), kmc(3, 0);
    for (int i = 0; i <= 2; ++i){
        for (int j = 0; j <= 2; ++j){
            kma[i] += khat[j] * Ma[j][i];
            kmb[i] += khat[j] * Mb[j][i];
            kmc[i] += khat[j] * Mc[j][i];
        }
    }
    polA = polarAngle(kma);
    azimA = azimuthalAngle(kma);
    polB = polarAngle(kmb);
    azimB = azimuthalAngle(kmb);
    polC = polarAngle(kmc);
    azimC = azimuthalAngle(kmc);
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            YlmStar[0] = std::conj(Y_lm(l, m, polangle[0], azimangle[0]));
            coeffTerm += YlmStar[0];
            YlmStar[1] = std::conj(Y_lm(l, m, polA, azimA));
            coeffTerm += pow((-1), (nvec[0]-nvec[1]))*YlmStar[1];
            YlmStar[2] = std::conj(Y_lm(l, m, polB, azimB));
            coeffTerm += pow((-1), (nvec[1]-nvec[2]))*YlmStar[2];
            YlmStar[3] = std::conj(Y_lm(l, m, polC, azimC));
            coeffTerm += pow((-1), (nvec[2]-nvec[0]))*YlmStar[3];
            coeffVec.push_back(std::pow(imu, l) * coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE7(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 1;
    std::vector< std::complex<double> > YlmStar(N);
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Mb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Mb[n]);
        azimangle[n] = azimuthalAngle(k_Mb[n]);
    }
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            for (int j = 0; j < N; ++j){
                YlmStar[j] = std::conj(Y_lm(l, m, polangle[j], azimangle[j]));
                coeffTerm += pow((-1),j*(nvec[0] + nvec[1])) * YlmStar[j];
            }
            coeffVec.push_back(std::pow(imu, l) * coeffTerm);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE8(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    int N = 2;
    std::vector< std::complex<double> > YlmStar(4);
    std::complex<double> coeffTerm = 0;
    std::vector<double> polangle(4), azimangle(4);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Ma(N, std::vector<double>(3, 0)), k_Mb(N, std::vector<double>(3, 0)), k_MaMb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = 1; Mb[2][2] = 1;
    std::vector<double> kma(3, 0), kmb(3, 0), kmamb(3, 0);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            kma[i] += khat[j] * Ma[j][i];
        }
    }
    polangle[1] = polarAngle(kma);
    azimangle[1] = azimuthalAngle(kma);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            kmb[i] += khat[j] * Mb[j][i];
        }
    }
    polangle[2] = polarAngle(kmb);
    azimangle[2] = azimuthalAngle(kmb);
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            kmamb[i] += khat[j] * Mb[j][i] * Ma[j][i];
        }
    }
    polangle[3] = polarAngle(kmamb);
    azimangle[3] = azimuthalAngle(kmamb);
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            for (int j = 0; j <= 3 ; ++j){
                YlmStar[j] = Y_lm(l, m, polangle[j], azimangle[j]);
            }
            coeffVec.push_back(std::pow(imu, l) * (YlmStar[0] )
                + pow((-1),(nvec[0] + nvec[1])) * YlmStar[1]
                + pow((-1),nvec[2]) * YlmStar[2]
                + pow((-1),(nvec[0]+nvec[1]+nvec[2])) * YlmStar[3]);
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE9(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    std::vector< std::complex<double> > YlmStar(2);
    std::complex<double> coeffTerm = 0;
    int N = 1;
    std::vector<double> polangle(N), azimangle(N);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Mb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > M(3, std::vector<double>(3, 0));
    M[0][0] = 1; M[1][1] = -1; M[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Mb[n]);
        azimangle[n] = azimuthalAngle(k_Mb[n]);
    }
    for (int l = 0; l <= lMax; ++l){
        for (int m = -l; m < l + 1; ++m){
            for (int j = 0; j < N; ++j){
                YlmStar[j] = std::conj(Y_lm(l, m, polangle[j], azimangle[j]));
                coeffTerm += pow((-1),j*(nvec[0] + nvec[1])) * YlmStar[j];
            }
            coeffVec.push_back((std::pow(imu, l) * coeffTerm));
        }
    }

    return coeffVec;
}

std::vector< std::complex<double> > CE10(std::vector<double> kvec, std::vector<int> nvec, std::vector<double> lengths, int lMax){
    std::vector< std::complex<double> >  coeffVec;
    std::vector<double> khat = hat(kvec);
    std::vector< std::complex<double> > YlmStar(2);
    int N = 1;
    std::vector<double> polangle(4), azimangle(4);
    polangle[0] = polarAngle(khat);
    azimangle[0] = azimuthalAngle(khat);
    std::vector< std::vector<double> > k_Ma(N, std::vector<double>(3, 0)), k_Mb(N, std::vector<double>(3, 0)), k_MaMb(N, std::vector<double>(3, 0));
    std::vector< std::vector<double> > Ma(3, std::vector<double>(3, 0));
    Ma[0][0] = 1; Ma[1][1] = -1; Ma[2][2] = 1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1; Mb[1][1] = -1; Mb[2][2] = 1;
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Ma[n][i] += khat[j] * matPow(Ma,n)[j][i];
            }
        }
        polangle[n] = polarAngle(k_Ma[n]);
        azimangle[n] = azimuthalAngle(k_Ma[n]);
    }
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_Mb[n][i] += khat[j] * matPow(Mb,n)[j][i];
            }
        }
    }
    polangle[2] = polarAngle(k_Mb[0]);
    azimangle[2] = azimuthalAngle(k_Mb[0]);
    for (int n = 0; n < N; ++n){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                k_MaMb[n][i] += khat[j] * matPow(Ma,n)[j][i] * matPow(Mb,n)[j][i];
            }
        }
    }
    polangle[3] = polarAngle(k_MaMb[0]);
    azimangle[3] = azimuthalAngle(k_MaMb[0]);
    for (int l = 2; l <= lMax; ++l){
        for(int m = -l; m < l + 1; ++m){
            for (int j = 0; j <= 3 ; ++j){
                YlmStar[j] = std::conj(Y_lm(l, m, polangle[j], azimangle[j]));
            }
            coeffVec.push_back(std::pow(imu, l) * (YlmStar[0] )
                + pow((-1),(nvec[0] + nvec[1])) * YlmStar[1]
                + pow((-1),nvec[2]) * YlmStar[2]
                + pow((-1),(nvec[0]+nvec[1]+nvec[2])) * YlmStar[3]);
        }
    }

    return coeffVec;
}

/**==============================__Normalization_Constants__===================================**/
/// some notation:
/// NEI -> normalization constant for Ith fundamental domain's eigenmode basis

/// inputs:
/// kvec: 3D wave-vector; must match domIndx of kvec to domIndx of CEI; must match domIndx to NEI as well
/// coeffs: coefficients of the eigenmode's basis; outputs of CEI functions


std::complex<double> NE1(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    for (int a = 0; a < 3; ++a){
        N += coeffs[0] * conj(coeffs[0]);// - imu * coeffs[0] * conj(coeffs[0]) * kvec[a];
    }//a
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE2(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 1;
    std::vector< std::vector<double> > Mb(3, std::vector<double>(3, 0));
    Mb[0][0] = -1.0; Mb[1][1] = -1.0; Mb[2][2] = 1.0;
    std::vector< std::vector<double> > I(3, std::vector<double>(3, 0));
    I[0][0] = 1.0; I[1][1] = 1.0; I[2][2] = 1.0;
    N += coeffs[0] * conj(coeffs[0]);
    N += coeffs[1] * conj(coeffs[1]);
    for (int a = 0; a < 3; ++a){
        for (int b = 0; b < 3; ++b){
            N += coeffs[0] * conj(coeffs[1]) + imu * coeffs[0] * conj(coeffs[1]) * kvec[a] * (Mb[a][b] - I[a][b]);
            N += coeffs[1] * conj(coeffs[0]) + imu * coeffs[1] * conj(coeffs[0]) * kvec[a] * (I[a][b] - Mb[a][b]);
        }//b
    }//a
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE3(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 3;
    std::vector< std::vector< std::vector<double> > > Mpows (q_max+1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //Mb
    Mpows[1][0][1] = 1; Mpows[1][1][1] = -1; Mpows[1][2][2] = 1; //Mb
    Mpows[2] = matPow(Mpows[1], 2); //Mb^2
    Mpows[3] = matPow(Mpows[1], 3); //Mb^3
    for (int q = 0; q <= q_max; ++q){
        for (int qp = 0; qp <= q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);
                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE4(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 2;
    std::vector< std::vector< std::vector<double> > > Mpows (q_max+1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = -0.5; Mpows[1][0][1] = -sqrt(3)/2; Mpows[1][1][0] = sqrt(3)/2; Mpows[1][1][1] = -0.5; Mpows[1][2][2] = 1; //Mb
    Mpows[2] = matPow(Mpows[1], 2); //Mb^2
    for (int q = 0; q <= q_max; ++q){
        for (int qp = 0; qp <= q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE5(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 5;
    std::vector< std::vector< std::vector<double> > > Mpows (q_max+1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = 0.5; Mpows[1][0][1] = -sqrt(3)/2; Mpows[1][1][0] = sqrt(3)/2; Mpows[1][1][1] = 0.5; Mpows[1][2][2] = 1; //Mb
    Mpows[2] = matPow(Mpows[1], 2); //Mb^2
    Mpows[3] = matPow(Mpows[1], 3); //Mb^3
    Mpows[4] = matPow(Mpows[1], 4); //Mb^4
    Mpows[5] = matPow(Mpows[1], 5); //Mb^5
    for (int q = 0; q <= q_max; ++q){
        for (int qp = 0; qp <= q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE6(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 4;
    std::vector< std::vector< std::vector<double> > > Mpows (4, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = 1; Mpows[1][1][1] = -1; Mpows[1][2][2] = -1; //Ma
    Mpows[2][0][0] = -1; Mpows[2][1][1] = 1; Mpows[2][2][2] = -1; //Mb
    Mpows[3][0][0] = -1; Mpows[3][1][1] = -1; Mpows[3][2][2] = 1; //Mc
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);
                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE7(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 2;
    std::vector< std::vector< std::vector<double> > > Mpows (q_max+1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = 1; Mpows[1][1][1] = -1; Mpows[1][2][2] = 1; //Mb
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE8(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 3;
    std::vector< std::vector< std::vector<double> > > Mpows (4, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = -1; Mpows[1][1][1] = 1; Mpows[1][2][2] = 1; //Ma
    Mpows[2][0][0] = 1; Mpows[2][1][1] = -1; Mpows[2][2][2] = 1; //Mb
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);
                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE9(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 2;
    std::vector< std::vector< std::vector<double> > > Mpows (q_max+1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = 1; Mpows[1][1][1] = -1; Mpows[1][2][2] = 1; //Mb
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

std::complex<double> NE10(std::vector<double> kvec, std::vector< std::complex<double> > coeffs){
    std::complex<double> N = 0.0;
    int q_max = 3;
    std::vector< std::vector< std::vector<double> > > Mpows (4, std::vector< std::vector<double> >(3, std::vector<double>(3, 0)));
    Mpows[0][0][1] = 1; Mpows[0][1][1] = 1; Mpows[0][2][2] = 1; //I
    Mpows[1][0][0] = -1; Mpows[1][1][1] = 1; Mpows[1][2][2] = 1; //Ma
    Mpows[2][0][0] = -1; Mpows[2][1][1] = -1; Mpows[2][2][2] = 1; //Mb
    for (int q = 0; q < q_max; ++q){
        for (int qp = 0; qp < q_max; ++qp){
            for (int a = 0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    N += coeffs[q] * conj(coeffs[qp]);
                    if (q != qp){
                        N += imu * coeffs[q] * conj(coeffs[qp]) * kvec[a] * (Mpows[q][a][b] - Mpows[qp][a][b]);
                    }
                }//b
            }//a
        }//qp
    }//q
    if (N == 0.0) return std::complex<double>(0.0, 0.0);
    else return 1.0 / sqrt(N);
}

/**==============================__3-Vectors_of_Integers__===================================**/

EigBase::retNs EigBase::domainNvec(){
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

/**==============================__Eigenmode_Basis_Computation__================================**/

/// computes summation needed for each eigenbasis term (terms w/in the summation over l) for a given fundamental domain
///     inputs: k->kvec; outputs of KEI's... coeffs->ouputs of CEI's... lMax as specified above
std::vector< std::complex<double> > basis(std::vector<double> k, double Rlss, std::vector<double> rT, std::vector< std::complex<double> > coeffs, std::vector< std::vector<double> > M){////, int lMax){
    std::vector< std::complex<double> > basis;
    double knorm = norm(k);
    std::complex<double> basis_sum;
    double rT_theta = polarAngle(rT);
    double rT_phi = azimuthalAngle(rT);
    for (int ell = 2; ell < lMax; ell++){
        for (int em = -ell; em < ell+1; em++){
            for (int i = 0; i < coeffs.size(); i++){
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        basis_sum += std::exp(imu*k[j]*M[j][k]*rT[k]);
                    }//k
                }//j
                    basis_sum += coeffs[i] * j_l(ell, knorm * Rlss) * Y_lm(ell, em, rT_theta, rT_phi);
                        if (std::isnan(std::real( basis_sum ))) std::cout << "jl = " << j_l(j, knorm) << std::endl;
                basis.push_back(basis_sum);
            }//i
        }//em
    }//ell
    return basis;
}

EigBase::retVals EigBase::eigBasis(double Rlss, std::vector<double> rT){
    /// terms needed to specify the eigenmode basis
    EigBase::retVals ret;
    std::vector< std::complex<double> > Cf; // CEI output
    std::complex<double> nrmC;              // NEI output
    std::vector<double> kvec(3, 0);         // KEI output
    std::vector< std::vector<int> > nvec;   // array of domainNvec().nsTr outputs
    std::complex<double> esum = 0;          // eigenmode basis summation term
    std::vector< std::complex<double> > xiterm(6, 0);

//    E1
    if (domIndx == 1){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE1(kvec, nvec[n], lengths);//, 1);
            nrmC = NE1(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum = basis(kvec, Rlss, rT, Cf, M, lMax)[l];
            }//l
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
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE2(kvec, nvec[n], lengths);
            nrmC = NE2(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E3
    if (domIndx == 3){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE3(kvec, nvec[n], lengths);
            nrmC = NE3(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E4
    if (domIndx == 4){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Khx(nvec[n], angles, lengths);
            Cf = CE4(kvec, nvec[n], lengths);
            nrmC = NE4(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E5
    if (domIndx == 5){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Khx(nvec[n], angles, lengths);
            Cf = CE5(kvec, nvec[n], lengths);
            nrmC = NE5(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E6
    if (domIndx == 6){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Krd(nvec[n], angles, lengths);
            Cf = CE6(kvec, nvec[n], lengths);
            nrmC = NE6(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
//    E7
    if (domIndx == 7){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE7(kvec, nvec[n], lengths);
            nrmC = NE7(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(esum);
    }
//    E8
    if (domIndx == 8){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE8(kvec, nvec[n], lengths);
            nrmC = NE8(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
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
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE9(kvec, nvec[n], lengths);
            nrmC = NE9(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(esum);
    }
//    E10
    if (domIndx == 10){
        esum = 0;
        nvec = domainNvec().nsTr;
        int nsz = nvec.size();
        for (int n = 0; n < nsz; n++){
            kvec = Kll(nvec[n], angles, lengths);
            Cf = CE10(kvec, nvec[n], lengths);//, condition);
            nrmC = NE10(kvec, Cf);
            for (int l = 0; l < lMax; l++){
                esum += basis(kvec, Rlss, rT, Cf, lMax)[l];
            }//l
            ret.coeffs.push_back(Cf);
            ret.k.push_back(kvec);
        }//n
        ret.basis.push_back(nrmC*esum);
    }
    return ret;
}

void EigBase::set_lMax(int lMaxIn){
    lMax = lMaxIn;
}

int EigBase::get_lMax(){
    return lMax;
}

void EigBase::set_nMax(int nMaxIn){
    nMax = nMaxIn;
}

int EigBase::get_nMax(){
    return nMax;
}

void EigBase::set_domIndx(int domIndxIn){
    domIndx = domIndxIn;
}

int EigBase::get_domIndx(){
    return domIndx;
}

void EigBase::set_angles(std::vector<double> anglesIn){
    angles[0] = anglesIn[0];
    angles[1] = anglesIn[1];
    angles[2] = anglesIn[2];
}

std::vector<double> EigBase::get_angles(){
    return angles;
}

void EigBase::set_lengths(std::vector<double> lengthsIn){
    lengths[0] = lengthsIn[0];
    lengths[1] = lengthsIn[1];
    lengths[2] = lengthsIn[2];
}

std::vector<double> EigBase::get_lengths(){
    return lengths;
}

