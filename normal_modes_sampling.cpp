/* This code calculates the normal modes from the coupled GOA model.
  To compile: g++ -o file file.cpp -llapack -lrefblas -lgfortran
  Dependence: LAPACK package
 (c) Xiang Sun 2015
*/
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

double Omega = 1; //primary mode frequency
double beta = 1;
double eta = 1;
const int n_omega = 200;
const double omega_max = 15;
const double DAcoupling = 0.1;
const double y_0 = 1.0; //shift of primary mode
const double d_omega = omega_max / n_omega;//SD evenly sampling rate
const double d_omega_eff = omega_max / n_omega; //for Jeffective SD
const double omega_c = 1; // the cutoff frequency for ohmic
const double pi =3.14159265358979324;
const double RT_2PI = sqrt(2*pi);
const double hbar = 1;
const int N = n_omega; //number of degrees of freedom
double J_omega_ohmic(double omega, double eta);
double GAUSS(int *seed);
double** Create_matrix(int row, int col);
extern "C" {
    void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double *A,
                const int &LDA, double *W, double *WORK, const int &LWORK,
                int &INFO);
}

int main (int argc, char *argv[]) {

    int i, j, a, b;
    int w;
    double omega;
    double t;
    int dim = n_omega;
    int col = dim, row = dim;
    double *eig_val = new double [row];
    double **matrix;
    matrix = Create_matrix(row,col);
    int lwork = 6*col, info;
    double *work = new double [lwork];
    double omega_nm[n_omega]; //normal mode frequencies
    double req_nm[n_omega]; //req of normal modes (acceptor shift)
    double c_nm[n_omega];//coupling strength of normal modes
    double **D_matrix;
    D_matrix = Create_matrix(n_omega, n_omega);
    double c_bath[n_omega];
    double **TT;
    TT = Create_matrix(n_omega, n_omega);
    int seed;
    seed = time(NULL);
    srand(seed);
    double R[n_omega];
    double V[n_omega];

    //secondary bath modes
    for (w = 1; w < n_omega; w++) {
        c_bath[w] = sqrt( 2.0 / pi * J_omega_ohmic(w*d_omega, eta) * d_omega * d_omega * w);
    }
    
    //********** BEGIN of Normal mode analysis ***********
    for (i = 0; i < n_omega; i++)
        for (j = 0; j <n_omega; j++) D_matrix[i][j] = 0;
    D_matrix[0][0] = Omega*Omega;
    for (w =1 ; w < n_omega ; w++) {
        D_matrix[0][0] += pow(c_bath[w]/(w*d_omega) ,2);
        D_matrix[0][w] = D_matrix[w][0] = c_bath[w];
        D_matrix[w][w] = pow(w*d_omega ,2);
    }
    
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            matrix[j][i] = D_matrix[i][j];
        }
    }

    dsyev_('V', 'L', col, matrix[0], col, eig_val, work, lwork, info);
    
    for (i = 0; i < dim; i++) omega_nm[i] = sqrt(eig_val[i]);
    
    for (i=0; i < dim; i++)
        for (j=0; j < dim; j++) TT[i][j] = matrix[i][j];

    for (i = 0; i < n_omega; i++) {
        req_nm[i] = 1 * TT[i][0];
        for (a = 1; a < n_omega; a++) req_nm[i] -= TT[i][a] * c_bath[a] / (a*d_omega * a*d_omega);
    }
    
    for (i = 0; i < n_omega; i++) {
        c_nm[i] = req_nm[i] * omega_nm[i] * omega_nm[i];
        req_nm[i] *= 2.0 * y_0;
    }
    
    //******** END of Normal mode analysis **************
    
    
    double sigma_x[n_omega];
    double sigma_p[n_omega];
    
    //Initial condition sampling (Gaussian)
    for (w = 1 ; w < n_omega; w++) {
        R[w] = GAUSS(&seed) * sigma_x[w];
        V[w] = GAUSS(&seed) * sigma_p[w];
    }
    
    return 0;
}

//************* subroutines *****************

double J_omega_ohmic(double omega, double etaa) {
    return etaa * omega * exp(-1 * omega / omega_c);
}

double** Create_matrix(int row, int col) {
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (int i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    return matrix;
}

// generate Normal distribution
double GAUSS(int *seed) {
    double A1 = 3.949846138;
    double A3 = 0.252408784;
    double A5 = 0.076542912;
    double A7 = 0.008355968;
    double A9 = 0.029899776;
    double SUM, R, R2, random, rannum;
    SUM = 0.0;
    
    for (int i=1; i<=12; i++) {
        random = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
        SUM += random;
    }
    R = (SUM - 6.0)/4.0;
    R2 = R*R;
    rannum = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R;
    return (rannum);
}



