/*  This code calculates Fermi Golden Rule rate with linear coupling (non-Condon)
    using Brownian oscillator model
    To compile: g++ -o EFGRL EFGRL.cpp -llapack -lrefblas -lgfortran
    (c) Xiang Sun 2015
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "r_1279.h"  // r1279 random number generator
using namespace std;

string nameapp = "b10e0.1_";
double beta = 10;//0.2;//1;//5;
const double eta = 0.1; //0.2;//1;//5;
const double DAcoupling = 0.1;

double Omega = 0.5; //primary mode freq
double y_0 = 1; //shift of primary mode

//for gaussian spectral density
const double sigma = 0.1;
const double omega_op = 1.0;

const int n_omega = 200;// 500 for low temp;//200;
const double d_omega = 0.05;// 0.02; for low temp////0.1;//0.002;for gaussian//0.1; for ohmic
//const double omega_max = 1;//20;//2.5 for gaussian// 20 for ohmic
const double d_omega_eff = 0.05;//0.02; for low temp //for effective SD sampling rate


const int LEN = 512;// 1024 for low temp//1024; //number of t choices 1024 for gaussian//512 for ohmic
const double DeltaT= 0.2;//0.2;//0.3; for gaussian//0.2 for ohmic //FFT time sampling interval
const double T0= -DeltaT*LEN/2;//-DeltaT*LEN/2+DeltaT/2;
const double pi=3.14159265358979324;
const double RT_2PI= sqrt(2*pi);
const double hbar = 1;

double S_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double S_omega_drude(double omega, double eta);//Drude spectral density
double S_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
double J_omega_ohmic(double omega, double eta);//bath Ohmic SD
double J_omega_ohmic_eff(double omega, double eta); //effective SD for Ohmic bath
void Integrand_exact(double omega, double t, double &re, double &im);
void Integrand_LSC(double omega, double t, double &re, double &im);
void Integrand_CAV(double omega, double t, double &re, double &im);
void Integrand_CD(double omega, double t, double &re, double &im);
void Integrand_W0(double omega, double t, double &re, double &im);
void Integrand_2cumu(double omega, double t, double &re, double &im);
void Integrand_Marcus(double omega, double t, double &re, double &im);
void Linear_exact(double omega, double t, double req, double &re, double &im);
void Linear_LSC(double omega, double t, double req, double &re, double &im);
void Linear_CAV(double omega, double t, double req, double &re, double &im);
void Linear_CD(double omega, double t, double req, double &re, double &im);
void Linear_W0(double omega, double t, double req, double &re, double &im);
void Linear_Marcus(double omega, double t, double req, double &re, double &im);
void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double Integrate(double *data, int n, double dx);
double Sum(double *data, int n);


extern "C" {
    void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double *A,
                const int &LDA, double *W, double *WORK, const int &LWORK,
                int &INFO);
}


int main (int argc, char *argv[]) {
    
    stringstream ss;
    string emptystr("");
    string filename;
    string idstr("");

    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT
	
	while (nn < LEN ) {
		mm++;
		nn *= 2;
	} //nn is the first 2^m that larger LEN
	
	double *corr1 = new double [nn];
	double *corr2 = new double [nn];
    
    double *corr1_orig = new double [nn]; //shifted origin to T0
    double *corr2_orig = new double [nn];
    
    double t;
    int i, j, a, b;
    double omega;
    int w; //count of omega
    double integ_re[n_omega];
    double integ_im[n_omega];
    
    long seed;
    seed = seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number 
    
    ofstream outfile;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    
    double Er=0;
    double a_parameter=0;
    double SD[n_omega];
    double J_eff[n_omega];
    
    
    //setting up spectral density
    for (w = 0; w < n_omega; w++) SD[w] = S_omega_ohmic(w*d_omega, eta); //Ohmic spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = S_omega_drude(w*d_omega, eta); //Drude spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = S_omega_gaussian(w*d_omega, eta, sigma, omega_op);
    for (w = 0; w < n_omega; w++) J_eff[w] = J_omega_ohmic_eff(w*d_omega_eff, eta);
    
    //outfile.open("S(omega).dat");
    //for (i=0; i< n_omega; i++) outfile << SD[i] << endl;
    //outfile.close();
    //outfile.clear();
    outfile.open("J_eff(omega).dat");
    for (i=0; i< n_omega; i++) outfile << J_eff[i] << endl;
    outfile.close();
    outfile.clear();
    
    /*
    double integrand[n_omega];
    for (w = 0; w < n_omega; w++) integrand[w] = SD[w] * w *d_omega;
    integrand[0]=0;
    Er = Integrate(integrand, n_omega, d_omega);
    
    for (w = 1; w < n_omega; w++) integrand[w] = SD[w] * w *d_omega * w*d_omega /tanh(beta*hbar* w * d_omega*0.5);
    integrand[0]=0;
    a_parameter = 0.5 * Integrate(integrand, n_omega, d_omega);
    cout << "Er = " << Er << endl;
    cout << "a_parameter = " << a_parameter << endl;
     */
    
    cout << "-------------EFGRL----------------" << endl;
    double shift = T0 / DeltaT;
    double N = nn;
    cout << "shift = " << shift << endl;
    double linear_accum_re;
    double linear_accum_im;
    double linear_re;
    double linear_im;
    double temp_re;
    double temp_im;
    double req_eff[n_omega];//req of effective SD
    
    for (w = 1; w < n_omega; w++) req_eff[w] = sqrt(8 * hbar * J_eff[w] / (pi * w * d_omega_eff*w * d_omega_eff*w));//eq min for each eff normal mode
    
    double c_bath[n_omega]; //secondary bath mode min shift coefficients
    for (w = 1; w < n_omega; w++) {
        c_bath[w] = sqrt( 2 / pi * J_omega_ohmic(w*d_omega, eta) * d_omega * d_omega * w);
    }
    
    
    //********** BEGIN of Normal mode analysis ***********
    
    //Dimension of matrix (Check this every time)
    int dim = n_omega;
    
    //-------- initiate LAPACK -----------
    int col = dim, row = dim;
    //Allocate memory for the eigenvlues
    double *eig_val = new double [row];
    //Allocate memory for the matrix
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    //Parameters for dsyev_ in LAPACK
    int lwork = 6*col, info;
    double *work = new double [lwork];
    //-------------------------------------
    
    double **TT_ns = new double* [col];
    TT_ns[0] = new double [col*row];
    for (i = 1; i < col; ++i)
        TT_ns[i] = TT_ns[i-1] + row;
    
    //transformation matrix: [normal mode]=[[TT_ns]]*[system-bath]
    //TT_ns * D_matrix * T_sn = diag, eigenvectors are row-vector of TT_ns
    double omega_nm[n_omega]; //normal mode frequencies
    double req_nm[n_omega]; //req of normal modes (acceptor shift)
    double c_nm[n_omega];//coupling strength of normal modes
    double S_array[n_omega];//Huang-Rhys factor for normal modes
   
    double D_matrix[n_omega][n_omega];// the Hessian matrix
    for (i=0; i< n_omega; i++) for (j=0; j<n_omega ;j++) D_matrix[i][j] = 0;
    D_matrix[0][0] = Omega*Omega;
    for (w =1 ; w < n_omega ; w++) {
        D_matrix[0][0] += pow(c_bath[w]/(w*d_omega) ,2);
        D_matrix[0][w] = D_matrix[w][0] = c_bath[w];
        D_matrix[w][w] = pow(w*d_omega ,2);
    }
    
    //cout << "Hession matrix D:" << endl;
    for (i=0; i < dim; i++) {
        for (j=0; j < dim; j++) {
            matrix[j][i] = D_matrix[i][j]; //switch i j to match with Fortran array memory index
            //cout << D_matrix[i][j] << " ";
        }
        //cout << endl;
    }

    //diagonalize matrix, the eigenvectors transpose is in result matrix => TT_ns.
    dsyev_('V', 'L', col, matrix[0], col, eig_val, work, lwork, info); //diagonalize matrix
    if (info != 0) cout << "Lapack failed. " << endl;
    
    for (i=0; i < dim; i++) omega_nm[i] = sqrt(eig_val[i]);
    
    outfile.open("normal_mode_freq.dat");
    for (i=0; i < dim; i++) outfile << omega_nm[i] << endl;
    outfile.close();
    outfile.clear();
    
    //cout << "eigen values = ";
    //for (i=0; i < dim; i++) cout << eig_val[i] <<"    ";
    //cout << endl;
    
    for (i=0; i < dim; i++)
        for (j=0; j < dim; j++) TT_ns[i][j] = matrix[i][j];

    /*
    double **Diag_matrix = new double* [col];
    Diag_matrix[0] = new double [col*row];
    for (i = 1; i < col; ++i)
        Diag_matrix[i] = Diag_matrix[i-1] + row;
     
    cout << "diagonalized Hessian matrix: " << endl;
    for (i=0; i < dim; i++) {
        for (j=0; j < dim; j++) {
            for (a=0; a < dim; a++)
                for (b=0; b < dim; b++) Diag_matrix[i][j] += TT_ns[i][a]*D_matrix[a][b]*TT_ns[j][b]; //TT_ns * D * T_sn
            cout << Diag_matrix[i][j] << "    " ;
        }
        cout << endl;
    }

     cout << endl;
     cout << "transformation matrix TT_ns (TT_ns * D * T_sn = diag, eigenvectors are row-vector of TT_ns): " << endl;
     for (i=0; i < dim; i++) {
         for (j=0; j < dim; j++) cout << TT_ns[i][j] << "    " ;
         cout << endl;
     }
    */
    
    // the coefficients of linear electronic coupling in normal modes (gamma[j]=TT_ns[j][0]*gamma_y), here gamma_y=1
    double gamma_nm[n_omega];
    for (i=0; i<n_omega; i++) gamma_nm[i] = TT_ns[i][0];
    
    //req of normal modes (acceptor shift)
    for (i=0; i<n_omega; i++) {
        req_nm[i] = 1 * TT_ns[i][0];
        for (a=1; a < n_omega; a++) req_nm[i] -= TT_ns[i][a] * c_bath[a] / (a*d_omega * a*d_omega);
    }
    
    outfile.open("Huang-Rhys.dat");
    for (i=0; i<n_omega; i++) {
        //tilde c_j coupling strength normal mode
        c_nm[i] = req_nm[i] * omega_nm[i] * omega_nm[i];
        req_nm[i] *= 2 * y_0;
        //discrete Huang-Rhys factor
        S_array[i] = omega_nm[i] * req_nm[i] * req_nm[i] /2;
        outfile << S_array[i] << endl;
    }
    outfile.close();
    outfile.clear();
    
    //******** END of Normal mode analysis **************
    
    
    
    // Note: here use min-to-min omega_{DA} as Fourier frequency
    
    
    
    //test cases: [1c] Eq FGR using continuous SD J(\omega)
    //Exact or LSC approximation
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] =0;
        integ_im[0] =0;
        for (w = 1; w < n_omega; w++) {
            omega = w * d_omega_eff;
            Integrand_LSC(omega, t, integ_re[w], integ_im[w]);
            integ_re[w] *= 4*J_eff[w]/pi/(omega*omega);
            integ_im[w] *= 4*J_eff[w]/pi/(omega*omega);
        }
        integral_re = Integrate(integ_re, n_omega, d_omega);
        integral_im = Integrate(integ_im, n_omega, d_omega);
        
        corr1[i] = exp(-1 * integral_re) * cos(integral_im);
        corr2[i] = -1 * exp(-1 * integral_re) * sin(integral_im);
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
   
    outfile.open((emptystr + "Exact_EFGR_cont_" + nameapp + ".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();

    
    
    
    //test cases: [1d] Eq FGR using discreitzed SD J(\omega)
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] = 0;
        integ_im[0] = 0;
        linear_accum_re = 0;
        linear_accum_im = 0;
        for (w = 0; w < n_omega; w++) {
            Integrand_exact(omega_nm[w], t, integ_re[w], integ_im[w]);
            integ_re[w] *= S_array[w];
            integ_im[w] *= S_array[w];
        }
        integral_re = Sum(integ_re, n_omega);
        integral_im = Sum(integ_im, n_omega);
        corr1[i] = exp(-1 * integral_re) * cos(integral_im);
        corr2[i] = -1 * exp(-1 * integral_re) * sin(integral_im);
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    outfile.open((emptystr + "Exact_EFGR_" + nameapp + ".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
   
    /*
    //[1c] exact quantum EFGR Linear coupling with continuous SD
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] = 0;
        integ_im[0] = 0;
        linear_accum_re = 0;
        linear_accum_im = 0;
        for (w = 1; w < n_omega; w++) {
            omega = w * d_omega_eff;
            Integrand_exact(omega, t, integ_re[w], integ_im[w]);
            integ_re[w] *= 4*J_eff[w]/pi/(omega*omega);
            integ_im[w] *= 4*J_eff[w]/pi/(omega*omega);
            
            Linear_exact(omega, t, req_eff[w], linear_re, linear_im);
            //Note: gamma_nm is for normal modes not for evenly sampled eff omega!!!
            //Thus, this result does not mean anything
            linear_accum_re += linear_re * gamma_nm[w] * gamma_nm[w];
            linear_accum_re += linear_im * gamma_nm[w] * gamma_nm[w];
        }
        integral_re = Integrate(integ_re, n_omega, d_omega);
        integral_im = Integrate(integ_im, n_omega, d_omega);
        temp_re = exp(-1 * integral_re) * cos(integral_im);
        temp_im = -1 * exp(-1 * integral_re) * sin(integral_im);
        corr1[i] = temp_re * linear_accum_re - temp_im * linear_accum_im;
        corr2[i] = temp_re * linear_accum_im + temp_im * linear_accum_re;
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    outfile.open("Exact_EFGRL_cont.dat");
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    */
    

    //[2d] exact quantum EFGR Linear coupling with discrete normal modes
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] = 0;
        integ_im[0] = 0;
        linear_accum_re = 0;
        linear_accum_im = 0;
        for (w = 0; w < n_omega; w++) {
            Integrand_exact(omega_nm[w], t, integ_re[w], integ_im[w]);
            integ_re[w] *= S_array[w];
            integ_im[w] *= S_array[w];
            
            Linear_exact(omega_nm[w], t, req_nm[w], linear_re, linear_im);
            linear_accum_re += linear_re * gamma_nm[w] * gamma_nm[w];
            linear_accum_re += linear_im * gamma_nm[w] * gamma_nm[w];
        }
        integral_re = Sum(integ_re, n_omega);
        integral_im = Sum(integ_im, n_omega);
        temp_re = exp(-1 * integral_re) * cos(integral_im);
        temp_im = -1 * exp(-1 * integral_re) * sin(integral_im);
        corr1[i] = temp_re * linear_accum_re - temp_im * linear_accum_im;
        corr2[i] = temp_re * linear_accum_im + temp_im * linear_accum_re;
    }
    
    outfile.open((emptystr + "Exact_EFGRL_t_re_" + nameapp + ".dat").c_str());
    for (i=0; i<LEN; i++) outfile << corr1[i] << endl;
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "Exact_EFGRL_t_im_" + nameapp + ".dat").c_str());
    for (i=0; i<LEN; i++) outfile << corr2[i] << endl;
    outfile.close();
    outfile.clear();
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    outfile.open((emptystr + "Exact_EFGRL_" + nameapp + ".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();

    //[3d] LSC EFGR Linear coupling with discrete normal modes
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] =0;
        integ_im[0] =0;
        linear_accum_re = 0;
        linear_accum_im = 0;
        for (w = 0; w < n_omega; w++) {
            Integrand_LSC(omega_nm[w], t, integ_re[w], integ_im[w]);
            integ_re[w] *= S_array[w];
            integ_im[w] *= S_array[w];
            
            Linear_LSC(omega_nm[w], t, req_eff[w], linear_re, linear_im);
            linear_accum_re += linear_re * gamma_nm[w] * gamma_nm[w];
            linear_accum_re += linear_im * gamma_nm[w] * gamma_nm[w];
        }
        integral_re = Sum(integ_re, n_omega);
        integral_im = Sum(integ_im, n_omega);
        temp_re = exp(-1 * integral_re) * cos(integral_im);
        temp_im = -1 * exp(-1 * integral_re) * sin(integral_im);
        corr1[i] = temp_re * linear_accum_re - temp_im * linear_accum_im;
        corr2[i] = temp_re * linear_accum_im + temp_im * linear_accum_re;
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT

    for(i=0; i<nn; i++) { //shift time origin 
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
 
    
    outfile.open((emptystr + "LSC_EFGRL_" + nameapp + ".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    

    
    
    
    
    
    cout << "Omega = " << Omega << endl;
    cout << "DeltaT = " << DeltaT << endl;
    cout << "N = " << LEN << endl;
    cout << "df = " << 1.0/LEN/DeltaT << endl;
    cout << "f_max = " << 0.5/DeltaT << endl;
    cout << "beta = " << beta << endl;
    cout << " eta = " << eta << endl;
    
    cout << "Done." << endl;



    // Deallocate memory
    delete [] eig_val;
    delete [] matrix[0];
    delete [] matrix;
    delete [] work;

    return 0;
}



/********* SUBROUTINE *************/


//spectral densities

double S_omega_ohmic(double omega, double eta) {
    // S_omega= sum_j omega_j * Req_j^2 / 2 hbar delta(\omega - \omega_j)
    return eta * omega * exp(-1 * omega);
}

double S_omega_drude(double omega, double eta) {
    return eta * omega /(1 + omega*omega);
}

double S_omega_gaussian(double omega, double eta, double sigma, double omega_op) {
    return   0.5 / hbar * eta * omega * exp(-(omega - omega_op)*(omega - omega_op)/(2*sigma*sigma))/RT_2PI/sigma;
}

double J_omega_ohmic(double omega, double eta) {
    //notice definition J(omega) is different from S(omega)
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return eta * omega * exp(-1 * omega);
}

double J_omega_ohmic_eff(double omega, double eta) {
    //(normal mode) effective SD for Ohmic bath DOF
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return eta * omega * pow(Omega,4) / ( pow(Omega*Omega - omega*omega, 2) + eta*eta*omega*omega);
}


//min-to-min energy as Fourier transform frequency

void Integrand_exact(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))/tanh(beta*hbar*omega*0.5);
    im = sin(omega*t);
    return;
}

void Integrand_LSC(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))/tanh(beta*hbar*omega*0.5);
    im = sin(omega*t);
    return;
}

void Integrand_W0(double omega, double t, double &re, double &im) {
    re = omega*omega*t*t*0.5/tanh(beta*hbar*omega*0.5);
    im = omega*t;
    return;
}

void Integrand_CAV(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = sin(omega*t);
    return;
}

void Integrand_CD(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = omega*t;
    return;
}

void Integrand_2cumu(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = omega*t;
    return;
}

void Integrand_Marcus(double omega, double t, double &re, double &im) {
    re = omega*t*t/(beta*hbar);
    im = omega*t;
    return;
}

void Linear_exact(double omega, double t, double req, double &re, double &im) {
    double cth = 1.0 / tanh(beta*hbar*omega*0.5);
    re = 0.5*hbar/omega*cth*cos(omega*t) + 0.25*req*req*((1-cos(omega*t))*(1-cos(omega*t)) - sin(omega*t)*sin(omega*t)*cth*cth));
    im = -0.5*hbar/omega*sin(omega*t) + 0.5*req*req*cth*(1-cos(omega*t))*sin(omega*t);
    return;
}

void Linear_LSC(double omega, double t, double req, double &re, double &im) {
    double cth = 1.0 / tanh(beta*hbar*omega*0.5);
    re = 0.5*hbar/omega*cth*cos(omega*t) - 0.25*req*req* sin(omega*t)*sin(omega*t)*cth*cth;
    im = -0.5*hbar/omega*sin(omega*t) + 0.25*req*req*cth*(1-cos(omega*t))*sin(omega*t);
    return;
}

void Linear_CAV(double omega, double t, double req, double &re, double &im) {
    re = (beta*cos(omega*t) - req*req*sin(omega*t)*sin(omega*t))/(beta*beta*omega*omega);
    im = 0;
    return;
}

void Linear_CD(double omega, double t, double req, double &re, double &im) {
    re = (beta*cos(omega*t) - req*req*sin(omega*t)*sin(omega*t))/(beta*beta*omega*omega);
    im = 0;
    return;
}

void Linear_W0(double omega, double t, double req, double &re, double &im) {
    double cth = 1.0 / tanh(beta*hbar*omega*0.5);
    re = cth * (0.5/omega - cth * 0.25*req*req*omega*omega*t*t);
    im = 0;
    return;
}


void Linear_Marcus(double omega, double t, double req, double &re, double &im) {
    re = 1.0/(beta*omega*omega) - req*req*t*t/(beta*beta);
    im = 0;
    return;
}


void FFT(int dir, int m, double *x, double *y)
{/*
    This code computes an in-place complex-to-complex FFT Written by Paul Bourke
    x and y are the real and imaginary arrays of N=2^m points.
    dir =  1 gives forward transform
    dir = -1 gives reverse transform
    Formula: forward
                N-1
                ---
            1   \           - i 2 pi k n / N
    X(n) =  -    >   x(k) e                       = forward transform
            1   /                                    n=0..N-1
                ---
                k=0
  
    Formula: reverse
                N-1
                ---
            1  \           i 2 pi k n / N
    X(n) = ---  >   x(k) e                  = reverse transform
            N  /                               n=0..N-1
                ---
                k=0
*/
    
    int n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    // Calculate the number of points
    n = 1;
    for (i=0;i<m;i++)
        n *= 2;
    
    // Do the bit reversal
    i2 = n >> 1; //i2 = (010 ...0)_2,second highest bit of n=(100 ...0)_2
    j = 0; //reversely bit accumulater from the second highest bit, i2.
    for (i=0;i<n-1;i++) {
        if (i < j) {
            tx = x[i]; //swap(i,j)
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        //to find the highest non-one bit, k, from the second highest bit
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k; //add 1 reversly
    }
    
    // Compute the Radix-2 FFT: Cooley-Tukey Algorithm
    c1 = -1.0; // c1+i*c2 = -1 = c^(i 2Pi/2) = W_2, def W_N^j = e^(i 2j*Pi/N)
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                //Butterfly calculation of x,y[i] and x,y[i1]:
                //t1+i*t2 =(u1+i*u2)(x[i1]+i*y[i2]) where u1+i*u2=W_N^j=e^(i 2j*Pi/N)
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            // i1+i*u2 *= c1+i*c2, or W_N
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        //c1+i*c2 = sqrt(c1+i*c2) eg. W_2 --> W_4 ...
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    // times STEPS*DeltaT forward FFT (time --> freq)
    /*if (dir == 1) {
     for (i=0; i<n; i++) {
     x[i] *= 1;//DeltaT;
     y[i] *= 1;//DeltaT;
     }
     }*/
    
    // Scaling for inverse transform
    
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    /*
     //for symmetrical FT,
     double sqn;
     sqn = sqrt(n);
     for (i=0;i<n;i++) {
     x[i] /= sqn;
     y[i] /= sqn;
     }
     */
    
    return;
}


double Integrate(double *data, int n, double dx){
    double I =0;
    I += (data[0]+data[n-1]);// divided by /2;
    for (int i=1; i< n-1; i++) {
        I += data[i];
    }
    I *= dx;
    return I;
}

double Sum(double *data, int n){
    double I = 0;
    for (int i=0; i< n; i++) {
        I += data[i];
    }
    return I;
}





