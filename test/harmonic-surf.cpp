/* This code calculates Fermi Golden Rule rate constant for charge transfer
   using multimode harmonic oscillator model, under different approximation levels
   (c) Xiang Sun 2014
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "r_1279.h"  // r1279 random number generator
using namespace std;

double beta = 4;//0.1;//4;

const double hbar = 1;
const double eta = 1;
const double DAcoupling = 0.1;
//for gaussian spectral density
const double sigma = 0.1;
const double omega_op = 1.0;
const double omega_max = 2.5; //20;
const double d_omega = 0.025;
const int n_omega = static_cast<int>(omega_max / d_omega);
const int MCN = 1000000; //Monte Carlo sampling number

const int N = n_omega; //number of degrees of freedom
const int LEN = 1024; //number of t choices in FFT
const int STEPS = 50; //for correlation function skip steps
const double DeltaT=0.3; //FFT time sampling interval
const double T0= -153.6; //FFT time sampling starting point
const double pi=3.14159265358979;
const double RT_2PI= sqrt(2*pi);

double DT=0.01; //MD time step
double DTSQ2 = DT * DT * 0.5;
double DT2 = DT/2.0;
double ABSDT= abs(DT);



void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double J_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double J_omega_drude(double omega, double eta);//another spectral density
double J_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
void Integrand_QM(double omega, double t, double &re, double &im);
double Integrate(double *data, int n, double dx);
double GAUSS(long *seed);
void MOVEA(double R[], double V[], double F[]);
void MOVEB(double V[], double F[]);
void force_avg(double R[], double F[], double omega[], double req[]);
void force_donor(double R[], double F[], double omega[], double req[]);
double DU(double R[], double omega[], double req[]);
double DUi(double R[], double omega[], double req[], int i);
int Job_finished(int &jobdone, int count, int total);



int main (int argc, char *argv[]) {

    int id; //jobid
    int i,j,k;
    stringstream ss;
    string emptystr("");
    string idstr="OLT_1e6_";
    
    double *SD = new double [n_omega];
    double omega;
    int w; //count of omega
    
    //setting up spectral density
    for (w = 0; w < n_omega; w++) SD[w] = J_omega_ohmic(w*d_omega, eta); //Ohmic spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_drude(w*d_omega, eta); //Drude spectral density
    //for (w = 0; w < n_omega; w++) SD[w] = J_omega_gaussian(w*d_omega, eta, sigma, omega_op);
    
    //outfile.open("SpectralDensity.dat");
    //for (i=0; i< n_omega; i++) outfile << SD[i] << endl;
    //outfile.close();
    //outfile.clear();
    
    //cout << "# of argument: " << argc-1 << endl;
    if (argc > 1) {
        ss << argv[1];
        ss >> id;
        idstr += ss.str();
        ss.clear();
    }
    cout << "----------------------------" << endl;
    cout << ">>> Job id # " << id << endl;
    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT

	while (nn < LEN ) {
		mm++;
		nn *= 2;
	} //nn is the first 2^m that larger than or equal to LEN

	double *corr1 = new double [nn];
	double *corr2 = new double [nn];

    double *corr1_orig = new double [nn]; //shifted origin to T0
    double *corr2_orig = new double [nn];

    double t;
    double shift = T0 / DeltaT;
    double NN = nn;

    int jobdone(0);

    double *integ_re = new double [n_omega];
    double *integ_im = new double [n_omega];

    long seed;
    seed = seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number

    ofstream outfile;

    double integral_re, integral_im;
    integral_re = integral_im = 0;

    double Er=0;
    double a_parameter=0;

    double *integrand = new double [n_omega];
    for (w = 0; w < n_omega; w++) integrand[w] = SD[w] * w *d_omega;
    integrand[0]=0;
    Er = Integrate(integrand, n_omega, d_omega);

    for (w = 1; w < n_omega; w++) integrand[w] = SD[w] * w*d_omega * w*d_omega /tanh(beta*hbar* w * d_omega*0.5);
    integrand[0]=0;
    a_parameter = 0.5 * Integrate(integrand, n_omega, d_omega);
    cout << "Er = " << Er << endl;
    cout << "a_parameter = " << a_parameter << endl;


    //min-to-min energy Fourier frequency
    cout << "DeltaT = " << DeltaT << endl;
    cout << "LEN = " << LEN << endl;
    cout << "df = " << 1.0/LEN/DeltaT << endl;
    cout << "f_max = " << 0.5/DeltaT << endl;
    cout << "beta = " << beta << endl;
    cout << " eta = " << eta << endl;

    // LSC approximation (numerical)
    double mc_re[LEN];
    double mc_im[LEN];
    double mc_re_D[LEN];
    double mc_im_D[LEN];
    double req[n_omega]; //acceptor equilibrium position
    double R[n_omega];  //position
    double V[n_omega];  //velocity
    double R_D[n_omega];  //position
    double V_D[n_omega];  //velocity
    double R0[n_omega];  //position
    double V0[n_omega];  //velocity
    double F[n_omega];  //force
    double F_D[n_omega];  //force
    double omega_array[n_omega];
    int LENMD = static_cast<int>(LEN*DeltaT/DT); //number of steps for MD
    int NMD;
    double *du_accum = new double [LENMD];
    double *du_accum_D = new double [LENMD];
    double sigma_x[n_omega];
    double sigma_p[n_omega];
    int MDlen;
    int tau;
    double integral_du[LEN];
    double integral_du_D[LEN];
    double *cuu = new double [LENMD/STEPS]; //noneq <du(0) du(t)>_avg
    double *cuu_D = new double [LENMD/STEPS]; //noneq <du(0) du(t)>_donor
    double du0;
    double t_array[LEN];  // FFT time variable t = T0, T0+DeltaT...
    double constant=1; //1/pow(pi*hbar,n_omega);

    req[0]=sigma_x[0]=sigma_p[0]=0;
    for (w = 1; w < n_omega; w++) {// (1/pi*h)^n Prod_omega tanh(beta*hbar*omega/2)
        omega_array[w] = w * d_omega;
        //constant *= tanh(beta*hbar*omega_array[w]*0.5);
        req[w] = sqrt(d_omega*SD[w]*2*hbar/omega_array[w]);
        sigma_x[w] = sqrt(hbar/(2*omega_array[w]*tanh(0.5*beta*hbar*omega_array[w])));
        sigma_p[w] = sigma_x[w] * omega_array[w];
    }
    
    for (i=0; i < LEN; i++) {//prepare FFT time variable
        t_array[i] = T0 + DeltaT * i;
        mc_re[i] = mc_im[i] = 0;
        mc_re_D[i] = mc_im_D[i] = 0;
        cuu[i] = cuu_D[i] = 0;
    }
        
    for (j=0; j < MCN; j++) { //Monte Carlo phase-space integration (R,P)
        for (w=1; w < n_omega; w++) {
            R0[w] = R[w] = R_D[w] = GAUSS(&seed) * sigma_x[w];//initial conf sampling
            V0[w] = V[w] = V_D[w] = GAUSS(&seed) * sigma_p[w];//initial momentum sampling
        }

        //Forward MD propagation
        if (t_array[LEN-1] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        NMD = static_cast<int>(abs(t_array[LEN-1]/DT));
        //cout << "NMD = " << NMD << endl;
        
        du0 = DU(R0, omega_array, req);
        //dynamics on average surface
        force_avg(R, F, omega_array, req);
        for (tau =0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_array, req);
            if (tau%STEPS == 0) {
                cuu[tau/STEPS] += du0 * du_accum[tau]; //noneq cuu avg
            }
            MOVEA(R, V, F);
            force_avg(R, F, omega_array, req);
            MOVEB(V, F);
        }
        
        //dynamics on donor surface
        force_donor(R_D, F_D, omega_array, req);
        
        for (tau =0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum_D[tau] = DU(R_D, omega_array, req);
            if (tau%STEPS == 0) {
                cuu_D[tau/STEPS] += du0 * du_accum_D[tau]; //noneq cuu donor
            }
            MOVEA(R_D, V_D, F_D);
            force_donor(R_D, F_D, omega_array, req);
            MOVEB(V_D, F_D);
        }
        
        for (i=0; i< LEN; i++) {
            if (t_array[i] >= 0) {
                MDlen = static_cast<int>(abs(t_array[i] / DT));
                integral_du[i] = Integrate(du_accum, MDlen, DT); //check ABSDT or DT?
                mc_re[i] += cos(integral_du[i]/hbar);
                mc_im[i] += sin(integral_du[i]/hbar);
                integral_du_D[i] = Integrate(du_accum_D, MDlen, DT); //check ABSDT or DT?
                mc_re_D[i] += cos(integral_du_D[i]/hbar);
                mc_im_D[i] += sin(integral_du_D[i]/hbar);
            }
        }
        
        //Backward MD propagation
        for (w=1;w < n_omega; w++) {
            R[w] = R0[w]; //restore initial condition
            V[w] = V0[w];
            R_D[w] = R0[w];
            V_D[w] = V0[w];
        }
        if (t_array[0] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        NMD = static_cast<int>(abs(t_array[0]/DT)); //>0
        //cout << "b NMD = " << NMD << endl;
        
        //dynamics on average surface
        force_avg(R, F, omega_array, req);
        for (tau =0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_array, req);
            MOVEA(R, V, F);
            force_avg(R, F, omega_array, req);
            MOVEB(V, F);
        }
        
        //dynamics on donor surface
        force_donor(R_D, F_D, omega_array, req);
        for (tau =0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum_D[tau] = DU(R_D, omega_array, req);
            MOVEA(R_D, V_D, F_D);
            force_donor(R_D, F_D, omega_array, req);
            MOVEB(V_D, F_D);
        }
        
        for (i=0; i< LEN; i++) {
            if (t_array[i] < 0) {
                MDlen = static_cast<int>(abs(t_array[i]/ DT)); // MDlen should >= 0
                integral_du[i] = Integrate(du_accum, MDlen, DT); //check ABSDT or DT?
                mc_re[i] += cos(integral_du[i]/hbar);
                mc_im[i] += sin(integral_du[i]/hbar);
                integral_du_D[i] = Integrate(du_accum_D, MDlen, DT); //check ABSDT or DT?
                mc_re_D[i] += cos(integral_du_D[i]/hbar);
                mc_im_D[i] += sin(integral_du_D[i]/hbar);
            }
        }
        Job_finished(jobdone, j, MCN);
    }

    //<du(0) du(t)_D/avg>_D for average- and donor-surface dynamics
    NMD = static_cast<int>(abs(t_array[LEN-1]/DT)/STEPS);
    cout << "Length of Cuu = " << NMD << endl;
    for (i=0;i<NMD;i++) {
        cuu[i] /= MCN;
        cuu_D[i] /= MCN;
    }
    
    outfile.open((emptystr+"Cuu_"+idstr+".dat").c_str());
    for (i=0; i< NMD; i++) outfile << cuu[i] << endl;
    outfile.close();
    outfile.clear();

    outfile.open((emptystr+"Cuu_D_"+idstr+".dat").c_str());
    for (i=0; i< NMD; i++) outfile << cuu_D[i] << endl;
    outfile.close();
    outfile.clear();


    //Analyze dynamics on average surface
    for (i=0; i< LEN; i++) { //Monte Carlo averaging
        mc_re[i] /= MCN;
        mc_im[i] /= MCN;
        corr1[i] = mc_re[i] * constant; // k(t) re
        corr2[i] = mc_im[i] * constant; // k(t) im
    }

    outfile.open((emptystr+"nlsc_t_re_"+idstr+".dat").c_str());
    for (i=0; i< LEN; i++) outfile << corr1[i] << endl;
    outfile.close();
    outfile.clear();

    outfile.open((emptystr+"nlsc_t_im_"+idstr+".dat").c_str());
    for (i=0; i< LEN; i++) outfile << corr2[i] << endl;
    outfile.close();
    outfile.clear();

    FFT(-1, mm, corr1, corr2);//notice its inverse FT

    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/NN) - corr2[i] * sin(-2*pi*i*shift/NN);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/NN) + corr1[i] * sin(-2*pi*i*shift/NN);
    }


    outfile.open((emptystr+"nlsc_re_"+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;

    outfile.close();
    outfile.clear();

    outfile.open((emptystr+"nlsc_im_"+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr2_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();

    //Analyze dynamics on donor surface
    for (i=0; i< LEN; i++) { //Monte Carlo averaging
        mc_re_D[i] /= MCN;
        mc_im_D[i] /= MCN;
        corr1[i] = mc_re_D[i] * constant; // k(t) re
        corr2[i] = mc_im_D[i] * constant; // k(t) im
    }

    outfile.open((emptystr+"nlsc_t_re_D_"+idstr+".dat").c_str());
    for (i=0; i< LEN; i++) outfile << corr1[i] << endl;
    outfile.close();
    outfile.clear();

    outfile.open((emptystr+"nlsc_t_im_D_"+idstr+".dat").c_str());
    for (i=0; i< LEN; i++) outfile << corr2[i] << endl;
    outfile.close();
    outfile.clear();

    FFT(-1, mm, corr1, corr2);//notice its inverse FT

    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/NN) - corr2[i] * sin(-2*pi*i*shift/NN);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/NN) + corr1[i] * sin(-2*pi*i*shift/NN);
    }


    outfile.open((emptystr+"nlsc_re_D_"+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;

    outfile.close();
    outfile.clear();

    outfile.open((emptystr+"nlsc_im_D_"+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr2_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();

    cout << "Done: job # " << idstr << endl;

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
    for (int i=1; i<N; i++) {
        xx = R[i] + DT*V[i] + DTSQ2*F[i];
        //pbc(xx, yy, zz);
        R[i] = xx;
        V[i] += DT2*F[i];
    }
    return;
}


void MOVEB(double V[], double F[]) {
	//always call MOVEB after call force** to update force F[]
    for (int i=1; i<N; i++) {
        V[i] += DT2*F[i];
    }
    return;
}

void force_avg(double R[], double F[], double omega[], double req[]) {
    //avg harmonic oscillator potential
    for (int i=1; i<N; i++) {
        F[i] = - omega[i] * omega[i] * (R[i]- req[i] * 0.5);
    }
    return;
}

void force_donor(double R[], double F[], double omega[], double req[]) {
    //donor harmonic oscillator potential
    for (int i=1; i<N; i++) {
        F[i] = - omega[i] * omega[i] * R[i];
    }
    return;
}

double DU(double R[], double omega[], double req[]) {
    double du=0;
    for (int i=1; i<N; i++) du += req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

double DUi(double R[], double omega[], double req[], int i) {
    double du=0;
    du = req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

int Job_finished(int &jobdone, int count, int total) {
	int tenpercent;
	tenpercent = static_cast<int> (10 * static_cast<double> (count)/ static_cast<double> (total) );
	if ( tenpercent > jobdone ) {
        jobdone = tenpercent;
        cout << "Job finished "<< jobdone <<"0%. " << endl;
	}

	return tenpercent;
}
