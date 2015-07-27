/* This code calculates equilibrium Fermi Golden Rule rate
 in Condon case, Brownian oscillator model
 compare with linearized semiclassical methods
 [special case with Jeff spectral density for ohmic bath modes]
 To compile: g++ -o EFGR_Jeff_num EFGR_Jeff_num.cpp
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

//*********** change parameter *********
const int MCN = 50000;//Monte Carlo sampling rate
double Omega = 0.2; //primary mode freq
double beta = 1;//0.2;//1;//5;
double eta = 1; //0.5;//1;//5;
const int n_omega = 1000; //normal cases 200
const double omega_max = 15;//15 or 20 for ohmic
const int LEN = 512;//1024;//512; //number of t choices
const double DeltaT = 0.2;//0.3;//0.2; //FFT time sampling interval
double DT=0.002; //MD time step
const double DAcoupling = 0.1;
//*********** **************** *********

const int INIT_omega = 1; //first omega index, =0 for normal modes, =1 for Jeff
const double y_0 = 1.0; //shift of primary mode
const double d_omega = omega_max / n_omega;//SD evenly sampling rate
const double d_omega_eff = omega_max / n_omega;//for effective SD sampling rate
const double omega_c = 1; //cutoff freq for ohmic
const int N = n_omega; //number of degrees of freedom
const double T0 = -DeltaT*(LEN*0.5);
const double hbar = 1;
const double pi =3.14159265358979324;
const double RT_2PI = sqrt(2*pi);
double DTSQ2 = DT * DT * 0.5;
double DT2 = DT/2.0;
double ABSDT = abs(DT);

//Declare Subroutines
void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double S_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double S_omega_drude(double omega, double eta);//another spectral density
double S_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
double J_omega_ohmic(double omega, double eta);//bath Ohmic SD
double J_omega_ohmic_eff(double omega, double eta); //effective SD for Ohmic bath
void Integrand_LSC(double omega, double t, double &re, double &im);
void Integrand_LSC_inh(double omega, double t, double &re, double &im);
void Integrand_CL_avg(double omega, double t, double &re, double &im);
void Integrand_CL_donor(double omega, double t, double &re, double &im);
void Integrand_2cumu(double omega, double t, double &re, double &im);
void Integrand_2cumu_inh(double omega, double t, double &re, double &im);
void Integrand_NE_exact(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
double Integrate(double *data, int n, double dx);
double Integrate_from(double *data, int sp, int n, double dx);
double Sum(double *data, int n);
double** Create_matrix(int row, int col);//new continuous 2D array in heap
double GAUSS(long *seed);
void MOVEA(double R[], double V[], double F[]);
void MOVEB(double V[], double F[]);
void force_avg(double R[], double F[], double omega[], double req[]);
double DU(double R[], double omega[], double req[]);
double DUi(double R[], double omega[], double req[], int i);
int Job_finished(int &jobdone, int count, int total, int startTime);


int main (int argc, char *argv[]) {
    
    int id;
    stringstream ss;
    string emptystr("");
    string nameapp("");
    string filename("");
    string idstr("");
    
    //cout << "# of argument: " << argc-1 << endl;
    if (argc > 1) {
        ss << argv[1];
        idstr += ss.str();
        ss >> id;
        ss.clear();
    }
    
    ss.str("");
    nameapp = "";
    ss << "Omega" << Omega << "_b" << beta << "e" << eta << "_" << MCN << "_";
    nameapp = ss.str();
    
    cout << ">>> Start Job id # " << id << " of num EFGR in Condon case." << endl;
    
    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT
    
    while (nn < LEN) {
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
    
    ofstream outfile;
    ofstream outfile1;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    
    double J_eff[n_omega]; // Jeff  for Omega << omega_c

    double shift = T0 / DeltaT;
    double c_eff[n_omega];//c from "effective" SD
    double req_eff[n_omega];//req from "effective" SD
    double Er_bath(0);//reorganization energy for Ohmic bath
    double Er_eff(0); //reorganization energy for effective SD
    double a_parameter_eff(0);
    
    int jobdone(0);
    int startTime;
    startTime = time(NULL);
    
    long seed;
    seed = seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number
    
    
    //setting up spectral density Jeff(omega)
    for (w = INIT_omega; w < n_omega; w++) J_eff[w] = J_omega_ohmic_eff(w*d_omega_eff, eta);//Jeff1

    //calculate Er and a parameter for Jeff
    Er_eff = 0;
    a_parameter_eff = 0;
    for (w = INIT_omega; w < n_omega; w++) {
        omega =  w * d_omega_eff;
        //eq min for each "eff" normal mode
        req_eff[w] = sqrt(8.0 * d_omega_eff * J_eff[w] / pi / pow(omega,3));
        //c_alpha for each "eff" normal mode
        c_eff[w] = sqrt( 2.0 / pi * J_eff[w] * d_omega_eff * omega);
        
        Er_eff += J_eff[w] / omega * d_omega_eff * 4.0/pi;
        //Er_eff_RRww += 0.5 * req_eff[w] * req_eff[w]  * w * d_omega_eff * w * d_omega_eff;
        
        a_parameter_eff += 0.25*pow(omega,3)*req_eff[w]*req_eff[w]/tanh(beta*hbar*omega*0.5);
    }
    
    cout << "Er_eff   = " << Er_eff << endl;
    cout << "a_parameter_eff = " << a_parameter_eff << endl;

    
    // numerical LSC approximation of EFGR using Monte Carlo
    int MDlen;
    int tau;
    int LENMD = static_cast<int>(LEN*DeltaT/DT); //number of steps for MD
    int NMD;
    
    double mc_re[LEN];
    double mc_im[LEN];
    double R0[n_omega];//wigner initial sampling
    double V0[n_omega];
    double R[n_omega]; //avg dynamics w/ wigner sampling (LSC)
    double V[n_omega];
    double F[n_omega];
    double *du_accum = new double [LENMD];
    double sigma_x[n_omega];//standard deviation of position
    double sigma_p[n_omega];//standard deviation of velocity
    double integral_du[LEN];
    double t_array[LEN];  // FFT time variable t = T0, T0+DeltaT...
    double omega_array[n_omega];
    
    //set standard deviation of initial sampling
    req_eff[0] = sigma_x[0] = sigma_p[0] = 0;
    for (w = INIT_omega; w < n_omega; w++) {
        omega_array[w] = w * d_omega_eff;
        sigma_x[w] = sqrt(hbar/(2*omega_array[w]*tanh(0.5*beta*hbar*omega_array[w])));
        sigma_p[w] = sigma_x[w] * omega_array[w];
    }
    
    for (i = 0; i < LEN; i++) {//prepare FFT time array
        t_array[i] = T0 + DeltaT * i;
        mc_re[i] = mc_im[i] = 0;
    }
    
    int NMD_forward = static_cast<int>(abs(t_array[LEN-1]/DT));
    int NMD_backward = static_cast<int>(abs(t_array[0]/DT));
    if(NMD_forward > NMD_backward) NMD = NMD_forward; //get max of NMD_forward and NMD_backward
    else NMD = NMD_backward;
    
    
    //Begin Monte Carlo importance sampling
    
    for (j = 0; j < MCN; j++) { //Monte Carlo phase-space integration (R,P)
        for (w = INIT_omega ; w < n_omega; w++) {
            R0[w] = R[w] = GAUSS(&seed) * sigma_x[w];//Wigner initial conf sampling
            V0[w] = V[w] = GAUSS(&seed) * sigma_p[w];//Wigner initial momentum sampling
        }
        
        // --->>> Forward MD propagation
        if (t_array[LEN-1] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        //NMD = NMD_forward; NMD now is max of NMD_forward and NMD_backward
        
        //dynamics on average surface + wigner sampling
        force_avg(R, F, omega_array, req_eff);
        for (tau = 0 ; tau < NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_array, req_eff);
            MOVEA(R, V, F);
            force_avg(R, F, omega_array, req_eff);
            MOVEB(V, F);
        }
        
        for (i = 0; i < LEN; i++) {
            if (t_array[i] >= 0) {
                MDlen = static_cast<int>(abs(t_array[i] / DT));
                integral_du[i] = Integrate(du_accum, MDlen, DT); //notice sign of DT
                mc_re[i] += cos(integral_du[i]/hbar);
                mc_im[i] += sin(integral_du[i]/hbar);
            }
        }
        
        // ---<<< Backward MD propagation
        for (w = INIT_omega; w < n_omega; w++) {
            R[w] = R0[w]; //restore initial condition
            V[w] = V0[w];
        }
        if (t_array[0] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        //NMD = NMD_backward; NMD now is max of NMD_forward and NMD_backward
        
        //dynamics on average surface + wigner sampling
        force_avg(R, F, omega_array, req_eff);
        for (tau = 0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_array, req_eff);
            MOVEA(R, V, F);
            force_avg(R, F, omega_array, req_eff);
            MOVEB(V, F);
        }
        
        for (i = 0; i < LEN; i++) {
            if (t_array[i] < 0) {
                MDlen = static_cast<int>(abs(t_array[i]/ DT)); // MDlen should >= 0
                integral_du[i] = Integrate(du_accum, MDlen, DT); //notice sign of DT
                mc_re[i] += cos(integral_du[i]/hbar);
                mc_im[i] += sin(integral_du[i]/hbar);
            }
        }
        Job_finished(jobdone, j, MCN, startTime);
    }
    
    //average Monte Carlo LSC (dynamics on average surface + wigner sampling)
    for (i = 0; i < LEN; i++) { //Monte Carlo averaging
        mc_re[i] /= MCN;
        mc_im[i] /= MCN;
        corr1[i] = mc_re[i]; //  k(t) re
        corr2[i] = mc_im[i]; //  k(t) im
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i = 0; i < nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/nn) - corr2[i] * sin(-2*pi*i*shift/nn);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/nn) + corr1[i] * sin(-2*pi*i*shift/nn);
    }
    
    
    outfile.open((emptystr+"num_LSC_EFGR_Jeff_"+nameapp+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    
    
    /*
    //=============case: analytical Eq FGR using continuous SD J_eff(\omega)=============
    //[a] Exact or LSC approximation using Jeff
    outfile1.open("Integral_Jeff.dat");
    
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
        integral_re = Integrate_from(integ_re, 1, n_omega, d_omega_eff);
        integral_im = Integrate_from(integ_im, 1, n_omega, d_omega_eff);
        //integral_re = Integrate(integ_re, n_omega, d_omega_eff);
        //integral_im = Integrate(integ_im, n_omega, d_omega_eff);
        
        outfile1 << integral_re << "\t" << integral_im << endl;
        
        corr1[i] = exp(-1 * integral_re) * cos(integral_im);
        corr2[i] = -1 * exp(-1 * integral_re) * sin(integral_im);
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    outfile.open("Exact_EFGR_Jeff.dat");
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    outfile1.close();
    outfile1.clear();
    */
    
    
    //-------------- Summary ----------------
    
    cout << "----- Parameters ------- " << endl;
    cout << "Omega (primary) = " << Omega << endl;
    cout << "n_omega = " << n_omega << endl;
    cout << "omega_max = " << omega_max << endl;
    cout << "beta = " << beta << endl;
    cout << "eta  = " << eta << endl;
    cout << "NMC = " << MCN << endl;
    cout << "LEN = " << LEN << endl;
    cout << "DeltaT = " << DeltaT << endl;
    cout << "--------- END of EFGR in Condon case --------" << endl;
    

    return 0;
}



/********* SUBROUTINE *************/


//spectral densities

double S_omega_ohmic(double omega, double etaa) {
    return etaa * omega * exp(-1 * omega / omega_c);
}

double S_omega_drude(double omega, double etaa) {
    return etaa * omega /(1 + omega*omega);
}

double S_omega_gaussian(double omega, double etaa, double sigma, double omega_op) {
    return   0.5 / hbar * etaa * omega * exp(-(omega - omega_op)*(omega - omega_op)/(2*sigma*sigma))/RT_2PI/sigma;
}

double J_omega_ohmic(double omega, double etaa) {
    //notice definition J(omega) is different from S(omega)
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return etaa * omega * exp(-1 * omega / omega_c);
}

double J_omega_ohmic_eff(double omega, double etaa) {
    //(normal mode) effective SD for Ohmic bath DOF
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return etaa * omega * pow(Omega,4) / ( pow(Omega*Omega - omega*omega, 2) + etaa*etaa*omega*omega);
}

//min-to-min energy as Fourier transform frequency
void Integrand_LSC(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))/tanh(beta*hbar*omega/2);
    im = sin(omega*t);
    return;
}

void Integrand_LSC_inh(double omega, double t, double &re, double &im) {
    re = omega*omega*t*t/2/tanh(beta*hbar*omega/2);
    im = omega*t;
    return;
}

void Integrand_CL_avg(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = sin(omega*t);
    return;
}

void Integrand_CL_donor(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = omega*t;
    return;
}

void Integrand_2cumu(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = omega*t;
    return;
}

void Integrand_2cumu_inh(double omega, double t, double &re, double &im) {
    re = omega*t*t/(beta*hbar);
    im = omega*t;
    return;
}

//noneq FGR

void Integrand_NE_exact(double omega, double tp, double tau, double shift, double req, double &re, double &im) {//including Huang-Rhys factor S_j
    re = omega*req*req*0.5*(1-cos(omega*tau))/tanh(beta*hbar*omega/2);
    im = omega*req*req*0.5*sin(omega*tau) + omega*req*shift* (sin(omega*tp) - sin(omega*tp-omega*tau));
    return;
}


double Integrate(double *data, int n, double dx){
    double I = 0;
    I += (data[0]+data[n-1])/2;
    for (int i=1; i< n-1; i++) {
        I += data[i];
    }
    I *= dx;
    return I;
}

double Integrate_from(double *data, int sp, int n, double dx){
    double I(0);
    I += (data[sp]+data[n-1])/2;
    for (int i=sp+1; i< n-1; i++) {
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
  X(n) = ----  >   x(k) e                       = forward transform
  1   /                                    n=0..N-1
  ---
  k=0
  
  Formula: reverse
  N-1
  ---
  1  \           i 2 pi k n / N
  X(n) =  ---  >   x(k) e                  = reverse transform
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



double** Create_matrix(int row, int col) {
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (int i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    return matrix;
}

// generate Normal distribution of random variable(rannum) from a uniform distribution ran()
double GAUSS(long *seed) {
    double A1 = 3.949846138;
    double A3 = 0.252408784;
    double A5 = 0.076542912;
    double A7 = 0.008355968;
    double A9 = 0.029899776;
    double SUM, R, R2, random, rannum;
    SUM = 0.0;
    
    for (int i=1; i<=12; i++) {
        //random = ran2(seed);
        random = r1279();
        SUM += random;
    }
    R = (SUM - 6.0)/4.0;
    R2 = R*R;
    rannum = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R;
    return (rannum);
}

// SUROUTINE TO PERFORM VELOCITY VERLET ALGORITHM
void MOVEA(double R[], double V[], double F[]) {
    double xx;
    for (int i= INIT_omega ; i<N; i++) {
        xx = R[i] + DT*V[i] + DTSQ2*F[i];
        //pbc(xx, yy, zz);
        R[i] = xx;
        V[i] += DT2*F[i];
    }
    return;
}


void MOVEB(double V[], double F[]) {
    //always call MOVEB after call force** to update force F[]
    for (int i = INIT_omega ; i<N; i++) {
        V[i] += DT2*F[i];
    }
    return;
}

void force_avg(double R[], double F[], double omega[], double req[]) {
    //avg harmonic oscillator potential
    for (int i = INIT_omega; i < N; i++) {
        F[i] = - omega[i] * omega[i] * (R[i]- req[i] * 0.5);
    }
    return;
}

double DU(double R[], double omega[], double req[]) {
    double du=0;
    for (int i = INIT_omega; i < N; i++) du += req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

double DUi(double R[], double omega[], double req[], int i) {
    double du=0;
    du = req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

int Job_finished(int &jobdone, int count, int total, int startTime) {
    int tenpercent;
    int currentTime;
    tenpercent = static_cast<int> (10 * static_cast<double> (count)/ static_cast<double> (total) );
    if ( tenpercent > jobdone ) {
        jobdone = tenpercent;
        currentTime = time(NULL);
        cout << "Job finished "<< jobdone <<"0 %. Time elapsed " << currentTime - startTime << " sec." << endl;
    }
    return tenpercent;
}

