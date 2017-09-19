#include "rotAngles.hpp"
#include <iostream>

double rotAngles(int angchoice, double euler, std::vector<double> vec){
    double rodRot[3][3];
    std::vector<double> vhat(3, 0);
    double norm = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if ( norm == 0 ){
        return 0;
    }
    else if (vec[0] == 0 && vec[1] == 0){
        return 0;
    }
    else{
        for (int i = 0; i < 3; i++){
            vhat[i] = (1 / norm) * vec[i];
        }
        std::vector<int> zhat(3);
        zhat[0]=0; zhat[1]=0; zhat[2]=1;
        std::vector<double> v(3);
        v[0] = vhat[1] * zhat[2] - vhat[2] * zhat[1];
        v[1] = vhat[2] * zhat[0] - vhat[0] * zhat[2];
        v[2] = vhat[0] * zhat[1] - vhat[1] * zhat[0];

        double vx[3][3];
        vx[0][0] = 0; vx[0][1] = -v[2]; vx[0][2] = v[1];
        vx[1][0] = -vx[0][1]; vx[1][1] = vx[0][0]; vx[1][2] = -v[0];
        vx[2][0] = -vx[0][2]; vx[2][1] = vx[1][2]; vx[2][2] = vx[0][0];

        double s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); //2-norm(v)
        double c = vhat[2]; //khat dot zhat
        double frac = ( 1 - c ) / ( s * s );
        for (int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                rodRot[i][j] = 1 + vx[i][j] * frac * vx[i][j] * vx[i][j];
            }
        }

        double chi = atan(-rodRot[0][2] / rodRot[1][2]);
        double psi = acos(rodRot[2][2]);
        if (angchoice == 1){return chi - euler;}
        else if (angchoice == 2){return psi - euler;}
        else {return 0;}
    }
}

