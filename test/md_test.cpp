/* This code test MD code time reversibility.
 (c) Xiang Sun 2014
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "r_1279.h"  // r1279 random number generator
using namespace std;

double beta = 0.1;

const double hbar = 1;

const double eta = 0.1;
const double DAcoupling = 0.1;
//for gaussian spectral density
const double sigma = 0.1;
const double omega_op = 1.0;

const double omega_max = 20;
const double d_omega = 0.002;

const int N = 1; //number of degrees of freedom
const int LEN = 1024; //number of t choices
const double DeltaT=0.3; //FFT time sampling interval
const double T0= -153.6; //FFT time sampling starting point
const double pi=3.14159265358979;
const double RT_2PI= sqrt(2*pi);

double DT=0.05; //MD time step
double DTSQ2 = DT * DT * 0.5;
double DT2 = DT/2.0;
double ABSDT= abs(DT);

const int MCN = 50000000; //Monte Carlo sampling number


void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double J_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double J_omega_drude(double omega, double eta);//another spectral density
double J_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
void Integrand_QM(double omega, double t, double &re, double &im);
double Integrate(double *data, int n, double dx);
double GAUSS(long *seed);
void MOVEA(double R[], double V[], double F[]);
void MOVEB(double V[], double F[]);
void force_avg(double R[], double F[], double omega, double req);
double DU(double R[], double omega, double req);



int main (int argc, char *argv[]) {
    
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
    int i,j,k;
    double omega;
    int w; //count of omega
    int n_omega = static_cast<int>(omega_max / d_omega);
    double *integ_re = new double [n_omega];
    double *integ_im = new double [n_omega];
    
    long seed;
    seed 	= seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number
    GAUSS(&seed);
    
    ofstream outfile;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    
    double Er=0;
    double a_parameter=0;
    double *SD = new double [n_omega];
    
    //setting up spectral density
    for (w = 0; w < n_omega; w++) SD[w] = J_omega_ohmic(w*d_omega, eta); //Ohmic spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_drude(w*d_omega, eta); //Drude spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_gaussian(w*d_omega, eta, sigma, omega_op);
   
    double shift = T0 / DeltaT;
    double N = nn;
    
    // LSC approximation (numerical)
    double prod_omega_re;
    double prod_omega_im;
    double mc_re;
    double mc_im;
    double req;
    double R[2]; //1=N
    double V[2];  //initial phase space point
    double F[2];  //force
    double *du_accum = new double [static_cast<int>(LEN*DeltaT/DT)];
    double sigma_x;
    double sigma_p;
    int MDlen;
    int tau;
    double integral_du;
    double constant=1/pow(pi*hbar,n_omega);
    
    omega = 1;
    req = 0.1;
    t = T0 + (LEN-1) * DeltaT;
    t /= 400;
    double times;
    times = time(NULL);
    
    if (t < 0) DT = - ABSDT; //for MD propagation direction
    else DT = ABSDT;
    DT2= 0.5 * DT;
    MDlen = static_cast<int>(t/DT);
    
    sigma_x = sqrt(hbar/(2*omega*tanh(0.5*beta*hbar*omega)));
    sigma_p = sigma_x * omega;
    
    cout << " MC sampling rate =" << MCN << endl;
    cout << " omega= " << omega << endl;
    cout << " t = " << t<< endl;
    
    //for(i=0; i< 100000; i++) cout << GAUSS(&seed)*sigma_x << endl;
    
    
    for (i=0; i< 6;i++) {
    
    mc_re = mc_im = 0;
    
    for (j=0; j < MCN; j++) { //Monte Carlo phase-space integration (R,P)
        for (k=0;k<N;k++) {
            R[k] = GAUSS(&seed)* sigma_x;//initial conf sampling
            V[k] = GAUSS(&seed)* sigma_p;//initial momentum sampling
        }
        
        du_accum[0] = DU(R, omega, req);
        force_avg(R, F, omega, req);
        for (tau =0 ; tau< MDlen; tau++) {
            MOVEA(R, V, F);
            force_avg(R, F, omega, req);
            MOVEB(V, F);
            //cout << R[0] << "   " << V[0] << endl;
            du_accum[tau] = DU(R, omega, req); //record DU
        }
        integral_du = Integrate(du_accum, MDlen, ABSDT);  //check ABSDT or DT?
        mc_re += cos(integral_du/hbar);
        mc_im += sin(integral_du/hbar);

    }
    mc_re /= MCN;
    mc_im /= MCN;
    
    cout << " i = " << i << endl;
    cout << "  mc_re = " << mc_re << endl;
    cout << "  mc_im = " << mc_im << endl;
    
    cout << "    time = " << time(NULL) - times << " s." << endl;
    times = time(NULL);
    }
    
    
    
    return 0;
}



/********* SUBROUTINE *************/


//spectral densities

double J_omega_ohmic(double omega, double eta) {
    return eta * omega * exp(-1 * omega);
}

double J_omega_drude(double omega, double eta) {
    return eta * omega /(1 + omega*omega);
}

double J_omega_gaussian(double omega, double eta, double sigma, double omega_op) {
    return   0.5 / hbar * eta * omega * exp(-(omega - omega_op)*(omega - omega_op)/(2*sigma*sigma))/RT_2PI/sigma;
}


//min-to-min energy as Fourier transform frequency
void Integrand_QM(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))/tanh(beta*hbar*omega/2);
    im = sin(omega*t);
    return;
}


double Integrate(double *data, int n, double dx){
    double I =0;
    I += (data[0]+data[n-1])/2;
    for (int i=1; i< n-1; i++) {
        I += data[i];
    }
    I *= dx;
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
    
    // Scaling for inverse transform
    
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    return;
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
    for (int i=0; i<N; i++) {
        xx = R[i] + DT*V[i] + DTSQ2*F[i];
        //pbc(xx, yy, zz);
        R[i] = xx;
        V[i] += DT2*F[i];
    }
    return;
}


void MOVEB(double V[], double F[]) {
	//always call MOVEB after call force** to update force F[]
    for (int i=0; i<N; i++) {
        V[i] += DT2*F[i];
    }
    return;
}

void force_avg(double R[], double F[], double omega, double req) {
    //avg harmonic oscillator potential
    for (int i=0; i<N; i++) {
        F[i] = - omega * omega * (R[i]-req*0.5);
    }
    return;
}

double DU(double R[], double omega, double req) {
    double du=0;
    for (int i=0; i<N; i++) du += req*omega*omega*R[i]-0.5*req*req*omega*omega;
    return du;
}

