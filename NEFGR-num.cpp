/* This code calculates Non-equilibrium Fermi Golden Rule rate 
   in Condon case using Brownian oscillator model
   compare with linearized semiclassical methods  
   To compile: g++ -o NEFGR-num NEFGR-num.cpp -llapack -lrefblas -lgfortran
   (c) Xiang Sun 2015  */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "r_1279.h"  // r1279 random number generator
using namespace std;

//*********** change parameter *********
const double beta = 5; //3;
const double eta = 0.5;  //3;
double omega_DA_fix = 0; //fixed omega_DA, with scan tp
double s = 1;    //Noneq. initial shift of parimary mode
double Omega = 0.5; //primary mode freq
const int MCN = 10000;//50000; //Monte Carlo sample rate
//*********** **************** *********

//double tp_fix = 5; //fixed t' for noneq FGR rate k(t',omega_DA) with scan omega_DA
double y_0 = 1; //shift of primary mode
const double DAcoupling = 0.1;
const double tp_max = 20; //scanning tp option, DeltaTau as step
const double Deltatp = 0.2; //time step for tp (0 ~ tp_max)
const double DeltaTau =0.002; //time slice for t' griding

const int n_omega = 200;
const int N = n_omega; //system degrees of freedom
const double omega_max = 15;//20;//2.5 for gaussian// 20 for ohmic
const double d_omega = omega_max / n_omega; //0.1;
const double d_omega_eff = d_omega;//0.1;

const int LEN = 512;//512;//1024; //number of t choices 1024 for gaussian//512 for ohmic
const double DeltaT=0.2;//0.2;//0.3; for gaussian//0.2 for ohmic //FFT time sampling interval
const double T0= -DeltaT*(LEN*0.5);//-DeltaT*LEN/2+DeltaT/2;

const double pi=3.14159265358979324;
const double RT_2PI= sqrt(2*pi);
const double hbar = 1;
//for gaussian spectral density
const double sigma = 0.1;
const double omega_op = 1.0;

//for numerical MD
double DT=0.0005; //MD time step 0.002
double DTSQ2 = DT * DT * 0.5;
double DT2 = DT/2.0;
double ABSDT= abs(DT);

//*****declare subroutines******
void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double S_omega_ohmic(double omega, double eta); //ohmic with decay spectral density
double S_omega_drude(double omega, double eta);//another spectral density
double S_omega_gaussian(double omega, double eta, double sigma, double omega_op);//gaussian spectral density
double J_omega_ohmic(double omega, double eta);//bath Ohmic SD
double J_omega_ohmic_eff(double omega, double eta); //effective SD for Ohmic bath
void Integrand_LSC(double omega, double t, double &re, double &im);
void Integrand_NE_exact(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
void Integrand_NE_CAV(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
void Integrand_NE_CD(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
void Integrand_NE_W0(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
void Integrand_NE_Marcus(double omega, double tp, double tau,  double shift, double req, double &re, double &im);
double Integrate(double *data, int n, double dx);
double Sum(double *data, int n);
double** Create_matrix(int row, int col);
double GAUSS(long *seed);
void MOVEA(double R[], double V[], double F[]);
void MOVEB(double V[], double F[]);
void force_avg(double R[], double F[], double omega[], double req[]);
void force_donor(double R[], double F[], double omega[], double req[]);
double DU(double R[], double omega[], double req[]);
double DUi(double R[], double omega[], double req[], int i);
int Job_finished(int &jobdone, int count, int total, int startTime);


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
    string nameapp("");
    int id;
    
    //cout << "# of argument: " << argc-1 << endl;
    if (argc > 1) {
        ss << argv[1];
        idstr += ss.str();
        ss >> id;
        ss.clear();
    }

    cout << ">>> Start Job id # " << id << " of num NEFGR in Condon case." << endl;
    
    int startTime;
    startTime = time(NULL);
    int jobdone(0);
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
    
    long seed;
    seed = seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number
    
    ofstream outfile; //output to file
    ofstream outfile1;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    
    double Er=0;
    double SD[n_omega];
    double J_eff[n_omega];
    
    int M; //time slice for tau = 0 --> tp
    int m; //index of tau
    double shift = T0 / DeltaT;
    double N = nn;
    //cout << "shift = " << shift << endl;
    double linear_accum_re;
    double linear_accum_im;
    double linear_re;
    double linear_im;
    double temp_re;
    double temp_im;
    double req_eff[n_omega];//req of effective SD
    
    //Dimension of matrix (Check this every time)
    int dim = n_omega;
    //-------- initiate LAPACK -----------
    int col = dim, row = dim;
    //Allocate memory for the eigenvlues
    double *eig_val = new double [row];
    //Allocate memory for the matrix
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (int i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    //Parameters for dsyev_ in LAPACK
    int lwork = 6*col, info;
    double *work = new double [lwork];
    //-------------------------------------
    
    double **TT_ns;//[n_omega][n_omega];
    TT_ns = Create_matrix(n_omega, n_omega);
    //transformation matrix: [normal mode]=[[TT_ns]]*[system-bath]
    //TT_ns * D_matrix * T_sn = diag, eigenvectors are row-vector of TT_ns
    double omega_nm[n_omega]; //normal mode frequencies
    double req_nm[n_omega]; //req of normal modes (acceptor shift)
    double c_nm[n_omega];//coupling strength of normal modes
    double S_array[n_omega];//Huang-Rhys factor for normal modes
    double **D_matrix;//[n_omega][n_omega];// the Hessian matrix
    D_matrix = Create_matrix(n_omega, n_omega);
    //double Diag_matrix[n_omega][n_omega]; //for testing diagonalization
    double gamma_nm[n_omega]; //linear coupling coefficient
    double shift_NE[n_omega]; //the s_j shifting for initial sampling
    
    double tp; //t' for noneq preparation
    M = static_cast<int> (tp_max/DeltaTau);
    int tp_index;
    double tau;
    double d_omega_DA = 2 * pi / LEN / DeltaT; //omega_DA griding size
    double omega_DA;
    double kre, kim;
    double sum(0);
    double kneq(0);
    
    ss.str("");
    nameapp = "";
    ss << "b" << beta;
    ss << "e" << eta << "_";
    nameapp = ss.str();
    
    //setting up spectral density
    //for (w = 0; w < n_omega; w++) J_eff[w] = J_omega_ohmic_eff(w*d_omega_eff, eta);
    //for (w = 0; w < n_omega; w++) SD[w] = S_omega_ohmic(w*d_omega, eta); //Ohmic spectral density
    //for (w = 1; w < n_omega; w++) req_eff[w] = sqrt(8 * hbar * J_eff[w] / (pi * w * d_omega_eff*w * d_omega_eff*w));//eq min for each eff normal mode
    
    double c_bath[n_omega]; //secondary bath mode min shift coefficients
    for (w = 1; w < n_omega; w++) {
        c_bath[w] = sqrt( 2 / pi * J_omega_ohmic(w*d_omega, eta) * d_omega * d_omega * w);
    }
    
    //********** BEGIN of Normal mode analysis ***********
    
    //construct Hessian matrix
    for (i=0; i< n_omega; i++) for (j=0; j<n_omega ;j++) D_matrix[i][j] = 0;
    D_matrix[0][0] = Omega*Omega;
    for (w =1 ; w < n_omega ; w++) {
        D_matrix[0][0] += pow(c_bath[w]/(w*d_omega) ,2);
        D_matrix[0][w] = D_matrix[w][0] = c_bath[w];
        D_matrix[w][w] = pow(w*d_omega ,2);
    }
    
    for (i=0; i < dim; i++) {
        for (j=0; j < dim; j++) {
            matrix[j][i] = D_matrix[i][j];
            //NOTE: switch i j to match with Fortran array memory index
        }
    }
    
    //diagonalize matrix, the eigenvectors transpose is in result matrix => TT_ns.
    dsyev_('V', 'L', col, matrix[0], col, eig_val, work, lwork, info); //diagonalize matrix
    if (info != 0) cout << "Lapack failed. " << endl;
    
    for (i=0; i < dim; i++) omega_nm[i] = sqrt(eig_val[i]);//normal mode freq

    for (i=0; i < dim; i++)
        for (j=0; j < dim; j++) TT_ns[i][j] = matrix[i][j];
    
    // the coefficients of linear electronic coupling in normal modes (gamma[j]=TT_ns[j][0]*gamma_y), here gamma_y=1
    for (i=0; i<n_omega; i++) gamma_nm[i] = TT_ns[i][0];
    //Noneq initial shift of each mode
    for (i=0; i<n_omega; i++) shift_NE[i] = s * gamma_nm[i];
    
    //req of normal modes (acceptor's potential energy min shift)
    for (i=0; i<n_omega; i++) {
        req_nm[i] = 1 * TT_ns[i][0];
        for (a=1; a < n_omega; a++) req_nm[i] -= TT_ns[i][a] * c_bath[a] / (a*d_omega * a*d_omega);
    }
    
    //outfile.open("Huang-Rhys.dat");
    for (i=0; i<n_omega; i++) {
        //tilde c_j coupling strength normal mode
        c_nm[i] = req_nm[i] * omega_nm[i] * omega_nm[i];
        req_nm[i] *= 2 * y_0;
        //discrete Huang-Rhys factor
        S_array[i] = omega_nm[i] * req_nm[i] * req_nm[i] /2;
        //outfile << S_array[i] << endl;
    }
    outfile.close();
    outfile.clear();
    
    //******** END of Normal mode analysis **************
    
    //cout << "Normal mode analysis done. " << endl;
    
    // Noneq LSC in Condon case using discreitzed J(\omega)
    //option: fix omega_DA, scan tp = 0 - tp_max
    
    omega_DA = omega_DA_fix; //fix omega_DA

    //ss.clear();
    ss.str("");
    ss << "s" << s << "w" << omega_DA_fix << "_" << MCN << "_";
    nameapp += ss.str();
    
    double R0[n_omega];//wigner initial sampling
    double V0[n_omega];
    double R[n_omega]; //donor dynamics starting with (R0, V0)
    double V[n_omega];
    double F[n_omega];
    double Rav[n_omega]; //average dynamics starting with (R(tp), - V(tp))
    double Vav[n_omega];
    double Fav[n_omega];
    int NMD_D; //steps of MD propagation on donor surface
    int NMD_AV; //steps of MD propagation on average surface
    
    int tp_index_max = static_cast<int> (tp_max / Deltatp);
    
    //allocate uncontineous 2D array with different column length
    double **C_re_accum = new double *[tp_index_max];
    double **C_im_accum = new double *[tp_index_max];
    for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
        tp = tp_index * Deltatp;
        M = static_cast<int> (tp/DeltaTau);
        C_re_accum[tp_index] = new double [M];
        C_im_accum[tp_index] = new double [M];
        for (m=0; m<M; m++) {
            C_re_accum[tp_index][m] = 0; //initialize C_re_accum[][]
            C_im_accum[tp_index][m] = 0; //initialize C_im_accum[][]
        }
    }
    double sigma_R[n_omega];// standard deviation of R0 and V0
    double sigma_V[n_omega];
    for (w = 0; w < n_omega; w++) {
        sigma_R[w] = sqrt(hbar/(2*omega_nm[w]*tanh(0.5*beta*hbar*omega_nm[w])));
        sigma_V[w] = sigma_R[w] * omega_nm[w];
    }


    NMD_D = static_cast<int> (Deltatp/ABSDT);  //MD steps for Deltatp
    NMD_AV = static_cast<int> (DeltaTau/ABSDT);//MD steps for DeltaTau
    double *du_accum = new double [NMD_AV];
    double sum_du(0);
    
    outfile.open((emptystr+"num_LSC_NEFGR_C_"+nameapp+idstr+".dat").c_str());
    if (!outfile.is_open()) cout << "Error when open file!" << endl;

    // *********  Start Monte Carlo phase-space integration (R0,P0) *********
    for (j = 0; j < MCN; j++) {
        
        //generate gaussian distribution of initial condition
        for (w = 0; w < n_omega; w++) {
            //Wigner initial conf sampling w/ noneq shift
            R0[w] = R[w] = GAUSS(&seed) * sigma_R[w] - shift_NE[w];
            //Wigner initial momentum sampling
            V0[w] = V[w] = GAUSS(&seed) * sigma_V[w];
        }
        
        //tp loop --------------------
        for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
            tp = tp_index * Deltatp;

            //propagate on donor surface for Deltatp
            force_donor(R, F, omega_nm, req_nm);
            for (a = 0 ; a < NMD_D; a++) {
                MOVEA(R, V, F);
                force_donor(R, F, omega_nm, req_nm);
                MOVEB(V, F);
            }
            
            for (w = 0; w < n_omega; w++) {
                Rav[w] = R[w];
                Vav[w] = -V[w];
                //notice propagation direction：sign of DT, and integral
                //here we use flip velocity and propagate forward, so DT > 0
            }

            // tau loop -------------------
            M = static_cast<int> (tp/DeltaTau);//update for each tp
            sum_du = 0;
            for (m = 0; m < M; m++) {
                tau = m * DeltaTau;
                
                //propagate on average surface for DeltaTau
                force_avg(Rav, Fav, omega_nm, req_nm);
                for (a = 0 ; a < NMD_AV; a++) {
                    MOVEA(Rav, Vav, Fav);
                    force_avg(Rav, Fav, omega_nm, req_nm);
                    MOVEB(Vav, Fav);
                    //record DU every DT: DU is sum of all frequencies
                    du_accum[a] = DU(Rav, omega_nm, req_nm);
                }
                sum_du += Sum(du_accum, NMD_AV) * ABSDT;//check sign
                C_re_accum[tp_index][m] += cos(sum_du / hbar);
                C_im_accum[tp_index][m] += sin(sum_du / hbar);

            } // end of tau loop
 
        } // end of tp loop
        
        Job_finished(jobdone, j, MCN, startTime);
        
    } // end of MC j loop
    
    // *********  END. Monte Carlo phase-space integration (R0,P0) *********
    
    for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
        tp = tp_index * Deltatp;
        M = static_cast<int> (tp/DeltaTau);
        for (m = 0; m < M; m++) {
            outfile << C_re_accum[tp_index][m]/MCN << endl;
            outfile << C_im_accum[tp_index][m]/MCN << endl;
        }
    }
    outfile.close();
    outfile.clear();
    
    
    /*
    //-------------- analytical NEFGR  ---------------
    outfile.open((emptystr+"ana_LSC_NEFGR_"+nameapp+idstr+".dat").c_str());
    outfile1.open((emptystr+"ana_LSC_P_NEFGR_"+nameapp+idstr+".dat").c_str());
    sum=0;
    kneq=0;
    for (tp = 0; tp < tp_max; tp += Deltatp) {
        kre = kim = 0;
        M = static_cast<int> (tp/DeltaTau);
        for (m=0; m<M; m++) {//tau index
            tau = m * DeltaTau;
            integ_re[0] = 0;
            integ_im[0] = 0;
            for (w = 0; w < n_omega; w++) {
                Integrand_NE_exact(omega_nm[w], tp, tau, shift_NE[w], req_nm[w], integ_re[w], integ_im[w]);
            }
            integral_re = Sum(integ_re, n_omega);// *DeltaTau;
            integral_im = Sum(integ_im, n_omega);// *DeltaTau;
            temp_re = exp(-1 * integral_re) * cos(integral_im);
            temp_im = exp(-1 * integral_re) * sin(-1 * integral_im);
            kre += temp_re * cos(omega_DA*tau) - temp_im * sin(omega_DA*tau);
            kim += temp_re * sin(omega_DA*tau) + temp_im * cos(omega_DA*tau);
        }
        kre *= DeltaTau;
        kim *= DeltaTau;
        kneq = kre*2*DAcoupling*DAcoupling;
        outfile << kneq << endl;
        sum += kneq * Deltatp;//probability of donor state
        //outfile1 << 1 - sum << endl; //1 - int dt' k(t')
        outfile1 << exp(-1*sum) << endl; //P = exp(- int dt' k(t'))
    }
    outfile.close();
    outfile.clear();
    
    outfile1.close();
    outfile1.clear();
    //-------------- END. analytical NEFGR  ---------------
    */

    
    //-------------- Summary ----------------
    
    cout << "   Parameters: " << endl;
    cout << "       beta = " << beta << endl;
    cout << "        eta = " << eta << endl;
    cout << "       fix omega_DA = " << omega_DA_fix << endl;
    cout << "       normal modes n_omega = " << n_omega << endl;
    cout << "       initial shift s = " << s << endl;
    cout << ">>> Job done. Total time elapsed " << time(NULL) - startTime << endl;
    return 0;
}



/********* SUBROUTINE *************/


//spectral densities

double S_omega_ohmic(double omega, double eta) {
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
void Integrand_LSC(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))/tanh(beta*hbar*omega/2);
    im = sin(omega*t);
    return;
}

void Integrand_NE_exact(double omega, double tp, double tau, double shift, double req, double &re, double &im) {//including Huang-Rhys factor S_j
    re = omega*req*req*0.5*(1-cos(omega*tau))/tanh(beta*hbar*omega/2);
    im = omega*req*req*0.5*sin(omega*tau) + omega*req*shift* (sin(omega*tp) + sin(omega*tau - omega*tp));
    return;
}


void Integrand_NE_CAV(double omega, double tp, double tau, double shift, double req, double &re, double &im) {
    re = req*req/beta*(1-cos(omega*tau));
    im = omega*req*req*0.5*sin(omega*tau) + omega*req*shift * (sin(omega*tp) + sin(omega*tau - omega*tp));
    return;
}

void Integrand_NE_CD(double omega, double tp, double tau, double shift, double req, double &re, double &im) {
    re = req*req/beta*(1-cos(omega*tau));
    im = omega*req*req*0.5 * omega*tau + omega*req*shift * (sin(omega*tp) + sin(omega*tau - omega*tp));
    return;
}

void Integrand_NE_W0(double omega, double tp, double tau, double shift, double req, double &re, double &im) {
    re = omega*req*req*0.5 / tanh(beta*hbar*omega/2) * omega*omega*tau*tau*0.5;
    im = omega*req*req*0.5 * omega*tau + omega*req*shift * cos(omega*tp) * omega*tau;
    return;
}

void Integrand_NE_Marcus(double omega, double tp, double tau, double shift, double req, double &re, double &im) {
    re = omega*req*req*0.5 * omega*tau*tau / beta;
    im = omega*req*req*0.5 * omega*tau + omega*req*shift * cos(omega*tp) * omega*tau;
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

double** Create_matrix(int row, int col) {//allocate continuous memory for 2D matrix
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
    for (int i = 0; i<N; i++) {
        xx = R[i] + DT*V[i] + DTSQ2*F[i];
        //pbc(xx, yy, zz);
        R[i] = xx;
        V[i] += DT2*F[i];
    }
    return;
}


void MOVEB(double V[], double F[]) {
    //always call MOVEB after call force** to update force F[]
    for (int i = 0; i<N; i++) {
        V[i] += DT2*F[i];
    }
    return;
}

void force_avg(double R[], double F[], double omega[], double req[]) {
    //avg harmonic oscillator potential
    for (int i = 0; i<N; i++) {
        F[i] = - omega[i] * omega[i] * (R[i]- req[i] * 0.5);
    }
    return;
}

void force_donor(double R[], double F[], double omega[], double req[]) {
    //donor harmonic oscillator potential
    for (int i = 0; i<N; i++) {
        F[i] = - omega[i] * omega[i] * R[i];
    }
    return;
}

double DU(double R[], double omega[], double req[]) {
    double du = 0;
    for (int i = 0; i<N; i++) du += req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

double DUi(double R[], double omega[], double req[], int i) {
    double du = 0;
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

















