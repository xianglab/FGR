/* This code calculates average and error bars of a serious input files
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

const int NFILE = 50;
const int BUNCH = 5;
int MIN = 1;
double beta = 1;//0.2;//1;//5;
const double eta = 1; //0.2;//1;//5;
string nameapp = "Omega1_b1e1_50000_";
double DeltaE_max = 10; // cutoff max of DeltaE/h in FFT for relative error

double Omega = 1;//0.5; //primary mode freq 0.2, 0.5, 1
double y_0 = 1;//sqrt(10);//10;//sqrt(10.0);//1.0; //shift of primary mode
const double omega_max = 15;//20;//15 or 20 for Jeff
const int n_omega = 1000;//10000;//200;//10000; //100;
const double d_omega = omega_max / n_omega;//0.1;
const int LEN = 512;//1024;//512; //number of t choices
const double DeltaT = 0.2;//0.3;//0.2; //FFT time sampling interval

const double T0= -DeltaT*LEN/2;//-DeltaT*LEN/2+DeltaT/2;
const double pi=3.14159265358979;
const double RT_2PI= sqrt(2*pi);
const double hbar = 1;
const double DAcoupling = 0.1;
const double df = 1.0/LEN/DeltaT; //frequency step after FFT
const double dx = df * 2 * pi;//angular frequency step after FFT
const int avg_cut = static_cast<int> (DeltaE_max/dx); //max index of freq-domain rate (after FFT)

void average_error(string filename, int length);
void average_lsc(string filename, int length, int min, int filenum, double *result);
void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double J_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double J_omega_drude(double omega, double eta);//another spectral density
double J_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
double J_omega_ohmic_eff(double omega, double etaa);//Jeff for GOA model
void Integrand_LSC(double omega, double t, double &re, double &im);
double Integrate(double *data, int n, double dx);
double Integrate_from(double *data, int sp, int n, double dx);

int main (int argc, char *argv[]) {

    int id; //jobid
    int i,j,k;
    stringstream ss;
    string filename;
    
    string idstr="";

    if (argc > 1) {
        ss << argv[1];
        idstr += ss.str();
        ss >> id;
        ss.clear();
    }

    double t_array[LEN];  // FFT time variable
    for (i=0; i < LEN; i++) {//prepare FFT time variable
        t_array[i] = T0 + DeltaT * i;
    }

    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT
	while (nn < LEN ) {
		mm++;
		nn *= 2;
	} //nn is the first 2^m that larger than or equal to LEN
    int NN=static_cast<int>(nn/2);

	double *corr1 = new double [nn];
	double *corr2 = new double [nn];
    
    double *corr1_orig = new double [nn]; //shifted origin to T0
    double *corr2_orig = new double [nn];
    
    double t;
    double omega;
    int w; //count of omega
    double *integ_re = new double [n_omega];
    double *integ_im = new double [n_omega];

    ofstream outfile;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    
    double Er=0;
    double a_parameter=0;
    //double *SD = new double [n_omega];
    double J_eff[n_omega];
    
    //setting up spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_ohmic(w*d_omega, eta); //Ohmic spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_drude(w*d_omega, eta); //Drude spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_gaussian(w*d_omega, eta, sigma, omega_op);
    for (w = 1; w < n_omega; w++) J_eff[w] = J_omega_ohmic_eff(w*d_omega, eta);//Jeff1
    
    double *integrand = new double [n_omega];
    /*
    for (w = 0; w < n_omega; w++) integrand[w] = SD[w] * w *d_omega;
    integrand[0]=0;
    Er = Integrate(integrand, n_omega, d_omega);
    
    for (w = 1; w < n_omega; w++) integrand[w] = SD[w] * w*d_omega * w*d_omega /tanh(beta*hbar* w * d_omega*0.5);
    integrand[0]=0;
    a_parameter = 0.5 * Integrate(integrand, n_omega, d_omega);
     */

    double shift = T0 / DeltaT;
    double N = nn;
    
    cout << "------------------------------------" << endl;
    cout << "DeltaT = " << DeltaT << endl;
    cout << "LEN = " << LEN << endl;
    cout << "df = " << df << endl;
    cout << "f_max = " << 0.5/DeltaT << endl;
    cout << "beta = " << beta << endl;
    cout << " eta = " << eta << endl;
    cout << "NFILE = " << NFILE << endl;
    cout << "BUNCH = " << BUNCH << endl;
    cout << "------------------------------------" << endl;


    
    //LSC approximation (ANALYTICAL EXACT RESULT)
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0; //zero padding
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] =0;
        integ_im[0] =0;
        for (w = 1; w < n_omega; w++) {
            omega = w * d_omega;
            Integrand_LSC(omega, t, integ_re[w], integ_im[w]);
            integ_re[w] *= 4*J_eff[w]/pi/(omega*omega);
            integ_im[w] *= 4*J_eff[w]/pi/(omega*omega);
        }
        integral_re = Integrate_from(integ_re, 1, n_omega, d_omega);
        integral_im = Integrate_from(integ_im, 1, n_omega, d_omega);
        
        corr1[i] = exp(-1 * integral_re) * cos(integral_im);
        corr2[i] = -1 * exp(-1 * integral_re) * sin(integral_im);
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    double *rate_ana = new double [nn/2];
    double *rate_num = new double [nn/2];
    for (i=0; i<nn/2; i++) rate_ana[i] = corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling;

    if (avg_cut > nn/2) cout << "Cut off radius too big" << endl;

    filename = "num_EFGR_Jeff_"+nameapp;
    outfile.open((filename+"Relative_ERR.dat").c_str());
    double relative_error;

    cout << "Average relative error:" << endl;
    
    if (NFILE % BUNCH) cout << "BUNCH cannot devide NFILE" << endl;
    else for (j= 1; j <= NFILE; j += BUNCH) {
        relative_error = 0;
        average_lsc(filename, NN, MIN, j, rate_num);
        for (i = 1; i <= avg_cut; i++) relative_error += abs((rate_num[i] - rate_ana[i])/rate_ana[i]);
        outfile << relative_error/avg_cut/j << endl;
        cout << relative_error/avg_cut/j << endl;
    }
    
    outfile.close();
    outfile.clear();



    return 0;
}


// ***************** SUBROUTINE *******************
void average_lsc(string filename, int length, int min, int filenum, double *result) {
    ifstream infile;
    int i,j;
    double *data1 = new double[length];
    double *data1_sq = new double[length];
    double *data1_av = new double[length];
    stringstream ss;
    string fnum;
    for (i=0; i< length; i++) data1[i]= data1_sq[i]= data1_av[i] =0;
    
    for (j= 0; j< filenum; j++) {
        ss.str("");
        ss << j+min;
        fnum = ss.str();
        infile.open((filename+fnum+".dat").c_str());
        for (i=0; i< length; i++) {
            infile >> data1[i];
            data1_av[i] += data1[i];
            data1_sq[i] += data1[i]*data1[i];
        }
        infile.close();
        infile.clear();
    }
    for (i=0; i<length;i++) {
        data1_av[i] /= filenum;
        data1_sq[i] /= filenum;
    }
    
    for (i=0; i< length; i++) result[i] = data1_av[i];

    delete [] data1;
    delete [] data1_av;
    delete [] data1_sq;
    
    return;
}


void average_error(string filename, int length) {
    ofstream outfile;
    ifstream infile;
    int i,j;
    double *data1 = new double[length];
    double *data1_sq = new double[length];
    double *data1_av = new double[length];
    stringstream ss;
    string fnum;
    for (i=0; i< length; i++) data1[i]= data1_sq[i]= data1_av[i] =0;

    for (j= 0; j< NFILE; j++) {
        ss.str("");
        ss << j+MIN;
        fnum = ss.str();
        infile.open((filename+fnum+".dat").c_str());
        for (i=0; i< length; i++) {
            infile >> data1[i];
            data1_av[i] += data1[i];
            data1_sq[i] += data1[i]*data1[i];
        }
        infile.close();
        infile.clear();
    }
    for (i=0; i<length;i++) {
        data1_av[i] /= NFILE;
        data1_sq[i] /= NFILE;
    }
    ss.str("");
    ss << NFILE;
    fnum = ss.str();

    outfile.open((filename+"AVG"+fnum+".dat").c_str());
    for (i=0; i< length; i++) outfile << data1_av[i] << endl;
    outfile.close();
    outfile.clear();

    outfile.open((filename+"ERR"+fnum+".dat").c_str());
    for (i=0; i< length; i++) outfile << sqrt((data1_sq[i] - data1_av[i]*data1_av[i])/NFILE) << endl;
    outfile.close();
    outfile.clear();

    delete [] data1;
    delete [] data1_av;
    delete [] data1_sq;

    return;
}

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

double J_omega_ohmic_eff(double omega, double etaa) {
    //(normal mode) effective SD for Ohmic bath DOF
    //J_omega = pi/2 * sum_a c_a^2 / omega_a delta(omega - omega_a)
    return etaa * omega * pow(Omega,4) *y_0*y_0 / ( pow(Omega*Omega - omega*omega, 2) + etaa*etaa*omega*omega);
}

void Integrand_LSC(double omega, double t, double &re, double &im) {
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

double Integrate_from(double *data, int sp, int n, double dx){
    double I(0);
    I += (data[sp]+data[n-1])/2;
    for (int i=sp+1; i< n-1; i++) {
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



