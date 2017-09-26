/** Author: Joshua Connel Phelan Osborne, 2017 CWRU Dpt Phys
    notes:
        1. whenever a change is made, run the command: cbp2make -in T3_correlations.cbp -out makefile (in the '<filepath>/T3_correlations' directory)
           before running the 'make' command again (same directory)
**/

#include <iostream>
#include <thread>
#include <mutex>
#include <future>
#include <vector>
#include <complex>

#include "correlation_matrices.hpp"
#include "eigenbasis_lin_T3.hpp"
#include "readToVec.hpp"

using namespace std;

void CxxWrite(vector< vector< complex<double> > > Cxx, string mat_fname, int domn, int i, int j, int m, int n, int dim1, int dim2){
    ofstream R_Cxx, I_Cxx; // Real and Imaginary parts of Cxx
    R_Cxx.open ("Outputs/R_" + mat_fname +"-"+ to_string(domn) +"-"+ to_string(i) +"-"+ to_string(j) +"-"+ to_string(m) +"-"+ to_string(n)); // naming convention for searching through files
    I_Cxx.open ("Outputs/I_" + mat_fname +"-"+ to_string(domn) +"-"+ to_string(i) +"-"+ to_string(j) +"-"+ to_string(m) +"-"+ to_string(n));
    for ( int row = 0; row < dim1; row++){ //modes or pixel numbers
        for ( int col = 0; col < dim2; col++){ //ditto
            R_Cxx << real(Cxx[row][col]) << " "; //write to file
            I_Cxx << imag(Cxx[row][col]) << " ";
        }
    }
    R_Cxx.close();
    I_Cxx.close();
}

int main()
{
///    read in data
    vector< vector<double> > transFuncData = readToVec("transFuncData.csv");      // contains: 1st row input k, 2nd-last rows transfer functions for each l
    vector< vector<double> > angleData = readToVec("angleData.csv");              // contains: rows of 3D arrays to vary the fundamental domain tilt angles; first element is x-y tilt: theta, second is y-z: phi, third z-x: psi
    vector< vector<double> > lengthData = readToVec("lengthData.csv");            // contains: edge length scales of the universal covering space
    vector< vector<double> > eulAngData = readToVec("eulerAngData.csv");          // contains: rows of 3D arrays of the three euler angles of the fundamental domain
    vector< vector<double> > positionData = readToVec("positionData.csv");        // contains: rows of 3D arrays of observer position w/in fundamental domain
//    vector< vector<double> > dataMask = readToVec("dataMask.csv");                // contains: masking information for comparison to Planck data
///    global parameter values
    int resolution = 3;     // resolution of HEALPix map
    int enmax = 2;          // maximum and/or min (-enmax) of the set of integers constraining the wave-vectors of a fundamental domain
    int ellmax = 5;         // maximum ell of the spherical harmonic expansion of the temperature functions over pixel space
    double Rlss = 4.23;        // radius of the last scattering surface in Gpc
    double A = 1;           // signal amplitude
    int nPix = 12 * pow(2, resolution)*pow(2, resolution); // number of pixels in the specified pixel map (=12*Nside^2)
///    declare classes and set private member variables that don't depend on loops
    eigenbasis_lin_T3 E;
    E.set_nMax(enmax);
    correlation_matrices M;
    M.set_res(resolution); M.set_lMax(ellmax); M.set_nMax(enmax);
///    loop over range of topological parameters
    eigenbasis_lin_T3::retVals eigB;
    vector< vector< complex<double> > > Cvv, Cpp, Cuu;
    for (int i = 0; i < angleData.size(); i++){
        Cvv.clear(); Cpp.clear(); Cuu.clear();
        E.set_angles(angleData[i]);
        for (int j = 0; j < lengthData.size(); j++){
            E.set_lengths(lengthData[j]);
            for (int m = 0; m < eulAngData.size(); m++){
                for (int n = 0; n < positionData.size(); n++){
                ///    compute the eigenbasis for each domain and correlation matrices
                    for (int domn = 1; domn <= 2; domn++){
                        cout << domn << endl;
//                        cout << lengthData[3][0] << endl;
                        E.set_domIndx(domn);
                        eigB = E.eigBasis(Rlss, positionData[n]);
                    ///    (l,m) x (l',m') mode to mode correlations; where (l,m) = v (\nu)
                        Cvv = M.mode_mode(A, Rlss, eigB.k, eigB.basis, transFuncData, eulAngData[m]);
                    ///    mu x mu' correlations; where mu indexes the transformed basis so Cuu has same dimensions as Cvv
//                        Cuu = M.maskCorMat(dataMask, Cvv);
                    ///    (pix no.) x (pix no.)' correlations; size depends on res
//                        Cpp = M.pix_pix(Cuu); // or M.pix_pix(Cvv) if Cuu is not necessary
                    ///    write matrices to file
//                        CxxWrite(Cpp, "Cpp", domn, i, j, m, n, Cpp.size(), Cpp[0].size());
                        CxxWrite(Cvv, "Cvv", domn, i, j, m, n, Cvv.size(), Cvv[0].size());
//                        CxxWrite(Cuu, "Cuu", domn, i, j, m, n, Cuu.size(), Cuu[0].size());
                    }//domn
                }//n
            }//m
        }//j
    }//i

    return 0;
}
