/* This code calculates Fermi Golden Rule rate with linear coupling (non-Condon)
    using multiple harmonic oscillator model, under different approximation levels
    numerically with Monte Carlo importance sampling
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
double beta = 1;//0.2;//1;//5;
double eta = 1; //0.2;//1;//5;
double Omega = 0.5; //primary mode freq
const int n_omega = 200; //normal cases 200
const double omega_max = 10;//15 or 20 for ohmic
double DT=0.002; //MD time step
//*********** **************** *********

double y_0 = 1; //shift of primary mode
const double DAcoupling = 0.1;
const double d_omega = omega_max / n_omega;// ~ 0.1 for ohmic
const double d_omega_eff = omega_max / n_omega; //for effective SD sampling rate

const int LEN = 512; //number of t choices 512 for normal case
const double DeltaT = 0.2;//0.2 (LEN=512) or 0.3 (LEN=2014) for ohmic //FFT time step
const double T0 = -DeltaT*LEN/2;
const double pi = 3.14159265358979;
const double RT_2PI= sqrt(2*pi);
const double hbar = 1;
double DTSQ2 = DT * DT * 0.5;
double DT2 = DT/2.0;
double ABSDT= abs(DT);
const int N = n_omega; //degrees of freedom

void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data
double S_omega_ohmic(double omega, double eta); //S(omega) spectral density
double J_omega_ohmic(double omega, double eta);//J(omega) bath Ohmic SD
double J_omega_ohmic_eff(double omega, double eta); //J effective SD for Ohmic bath
void Integrand_exact(double omega, double t, double &re, double &im);
void Integrand_LSC(double omega, double t, double &re, double &im);
void Integrand_CAV(double omega, double t, double &re, double &im);
void Integrand_CD(double omega, double t, double &re, double &im);
void Integrand_W0(double omega, double t, double &re, double &im);
void Integrand_Marcus(double omega, double t, double &re, double &im);
void Linear_exact(double omega, double t, double req, double &re, double &im);
void Linear_LSC(double omega, double t, double req, double &re, double &im);
void Linear_CAV(double omega, double t, double req, double &re, double &im);
void Linear_CD(double omega, double t, double req, double &re, double &im);
void Linear_W0(double omega, double t, double req, double &re, double &im);
void Linear_Marcus(double omega, double t, double req, double &re, double &im);
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
void Linear_num_LSC(double *omega, double *r0, double *rt, double *p0, double *gm, double &re, double &im);
int Job_finished(int &jobdone, int count, int total, int startTime);


extern "C" {
    void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double *A,
                const int &LDA, double *W, double *WORK, const int &LWORK,
                int &INFO);
}


int main (int argc, char *argv[]) {
    
    int id;
    stringstream ss;
    string emptystr("");
    string nameapp("");
    string filename;
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
    ss << "b" << beta;
    ss << "e" << eta;
    nameapp = ss.str();
    
    cout << ">>> Start Job id # " << id << " of num EFGRL in non-Condon case." << endl;

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
    
    ofstream outfile;
    
    double integral_re, integral_im;
    integral_re = integral_im = 0;
    double *integrand = new double [n_omega];
    double Er=0;
    double a_parameter=0;
    double SD[n_omega];
    double J_eff[n_omega];

    double shift = T0 / DeltaT;
    double N = nn;
    double linear_sum_re;
    double linear_sum_im;
    
    long seed;
    seed = seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */
    r1279(); //[0,1] random number
    
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
    double **D_matrix; // the Hessian matrix
    D_matrix = Create_matrix(n_omega, n_omega);
    double **TT_ns;
    TT_ns = Create_matrix(n_omega, n_omega);
    //transformation matrix: [normal mode]=[[TT_ns]]*[system-bath]
    //TT_ns * D_matrix * T_sn = diag, eigenvectors are row-vector of TT_ns
    double omega_nm[n_omega]; //normal mode frequencies
    double req_nm[n_omega]; //req of normal modes (acceptor shift)
    double c_nm[n_omega];//coupling strength of normal modes
    double S_array[n_omega];//Huang-Rhys factor for normal modes
    double c_bath[n_omega]; //secondary bath mode min shift coefficients
    double gamma_array[n_omega]; //linear coupling coefficient

    int jobdone(0);
    int startTime;
    startTime = time(NULL);
    
    
    //------- setting up spectral density ----------
    for (w = 0; w < n_omega; w++) SD[w] = S_omega_ohmic(w*d_omega, eta);
    for (w = 0; w < n_omega; w++) J_eff[w] = J_omega_ohmic_eff((w+1)*d_omega_eff, eta);
    for (w = 0; w < n_omega; w++) req_eff[w] = sqrt(8 * hbar * J_eff[w] / (pi * (w+1) * d_omega_eff* (w+1) * d_omega_eff*(w+1)));//eq min for each Jeff normal mode
    
    for (w = 1; w < n_omega; w++) {//ohmic bath mode coupling, essential to exact discrete modes
        c_bath[w] = sqrt( 2 / pi * J_omega_ohmic(w*d_omega, eta) * d_omega * d_omega * w);
    }
    
    
    //********** BEGIN of Normal mode analysis ***********
    
    for (i=0; i< n_omega; i++) for (j=0; j<n_omega ;j++) D_matrix[i][j] = 0;
    D_matrix[0][0] = Omega*Omega;
    for (w =1 ; w < n_omega ; w++) {
        D_matrix[0][0] += pow(c_bath[w]/(w*d_omega) ,2);
        D_matrix[0][w] = D_matrix[w][0] = c_bath[w];
        D_matrix[w][w] = pow(w*d_omega ,2);
    }

    for (i=0; i < dim; i++)
        for (j=0; j < dim; j++)
            matrix[j][i] = D_matrix[i][j]; //switch i j to match with Fortran array memory index
        
    //diagonalize matrix, the eigenvectors transpose is in result matrix => TT_ns.
    dsyev_('V', 'L', col, matrix[0], col, eig_val, work, lwork, info); //diagonalize matrix
    if (info != 0) cout << "Lapack failed. " << endl;
    
    for (i=0; i < dim; i++) omega_nm[i] = sqrt(eig_val[i]);
    
    for (i=0; i < dim; i++)
        for (j=0; j < dim; j++) TT_ns[i][j] = matrix[i][j];

    // the coefficients of linear electronic coupling in normal modes (gamma[j]=TT_ns[j][0]*gamma_y), here gamma_y=1
    for (i=0; i<n_omega; i++) gamma_array[i] = TT_ns[i][0];
    
    //req of normal modes (acceptor shift)
    for (i=0; i<n_omega; i++) {
        req_nm[i] = 1 * TT_ns[i][0];
        for (a=1; a < n_omega; a++) req_nm[i] -= TT_ns[i][a] * c_bath[a] / (a*d_omega * a*d_omega);
    }
    for (i=0; i<n_omega; i++) {
        //tilde c_j coupling strength normal mode
        c_nm[i] = req_nm[i] * omega_nm[i] * omega_nm[i];
        req_nm[i] *= 2 * y_0;
        //discrete Huang-Rhys factor
        S_array[i] = omega_nm[i] * req_nm[i] * req_nm[i] / 2.0;
    }
    //******** END of Normal mode analysis **************
    
    //exact reorganization energy Er for Marcus theory from normal modes
    Er = 0;
    a_parameter = 0;
    //for (i = 0; i < n_omega; i++) Er += 2.0 * c_nm[i] * c_nm[i] / (omega_nm[i] * omega_nm[i]);
    for (i = 0; i < n_omega; i++) Er += 0.5 * omega_nm[i] * omega_nm[i] * req_nm[i] * req_nm[i]; //S_array[i] * omega_nm[i];
    for (i = 0; i < n_omega; i++) a_parameter += 0.5 * S_array[i] * omega_nm[i] * omega_nm[i] /tanh(beta*hbar* omega_nm[i] *0.5);
    cout << "Er (normal mode) = " << Er << endl;
    cout << "a_parameter (normal mode) = "<< a_parameter << endl;
    

    
    /*
    //[A] analytical: exact quantum EFGR Linear coupling with discrete normal modes
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0;
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] = 0;
        integ_im[0] = 0;
        linear_sum_re = 0;
        linear_sum_im = 0;
        for (w = 0; w < n_omega; w++) {
            Integrand_exact(omega_nm[w], t, integ_re[w], integ_im[w]);
            integ_re[w] *= S_array[w];
            integ_im[w] *= S_array[w];
            
            Linear_exact(omega_nm[w], t, req_nm[w], linear_re, linear_im);
            linear_sum_re += linear_re * gamma_array[w] * gamma_array[w];
            linear_sum_im += linear_im * gamma_array[w] * gamma_array[w];
        }
        integral_re = Sum(integ_re, n_omega);
        integral_im = Sum(integ_im, n_omega);
        temp_re = exp(-1 * integral_re) * cos(integral_im);
        temp_im = -1 * exp(-1 * integral_re) * sin(integral_im);
        corr1[i] = temp_re * linear_sum_re - temp_im * linear_sum_im;
        corr2[i] = temp_re * linear_sum_im + temp_im * linear_sum_re;
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i=0; i<nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    outfile.open((emptystr + "Exact_EFGRL_" + nameapp + ".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    

    //[B] analytical: LSC approximation using normal modes
    for (i = 0; i < nn; i++) corr1[i] = corr2[i] = 0;
    for (i = 0; i < LEN; i++) {
        t = T0 + DeltaT * i;
        integ_re[0] = 0;
        integ_im[0] = 0;
        linear_sum_re = 0;
        linear_sum_im = 0;
        for (w = 0; w < n_omega; w++) {
            Integrand_LSC(omega_nm[w], t, integ_re[w], integ_im[w]);
            integ_re[w] *= S_array[w];
            integ_im[w] *= S_array[w];
            
            Linear_LSC(omega_nm[w], t, req_nm[w], linear_re, linear_im);
            linear_sum_re += linear_re * gamma_array[w] * gamma_array[w];
            linear_sum_im += linear_im * gamma_array[w] * gamma_array[w];
        }
        integral_re = Sum(integ_re, n_omega);
        integral_im = Sum(integ_im, n_omega);
        temp_re = exp(-1 * integral_re) * cos(integral_im);
        temp_im = -1 * exp(-1 * integral_re) * sin(integral_im);
        corr1[i] = temp_re * linear_sum_re - temp_im * linear_sum_im;
        corr2[i] = temp_re * linear_sum_im + temp_im * linear_sum_re;
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
     */
    
    // [C] numerical: LSC approximation of EFGRL
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
    double *linear_accum_re = new double [LENMD];
    double *linear_accum_im = new double [LENMD];
    double sigma_x[n_omega];//standard deviation of position
    double sigma_p[n_omega];//standard deviation of velocity
    double integral_du[LEN];
    double t_array[LEN];  // FFT time variable t = T0, T0+DeltaT...
    
    //set standard deviation of initial sampling
    for (w = 0; w < n_omega; w++) {// (1/pi*h)^n Prod_omega tanh(beta*hbar*omega/2)
        sigma_x[w] = sqrt(hbar/(2*omega_nm[w]*tanh(0.5*beta*hbar*omega_nm[w])));
        sigma_p[w] = sigma_x[w] * omega_nm[w];
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
        for (w = 0 ; w < n_omega; w++) {
            R0[w] = R[w] = GAUSS(&seed) * sigma_x[w];//Wigner initial conf sampling
            V0[w] = V[w] = GAUSS(&seed) * sigma_p[w];//Wigner initial momentum sampling
        }
        
        // --->>> Forward MD propagation
        if (t_array[LEN-1] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        //NMD = NMD_forward;

        //dynamics on average surface + wigner sampling
        force_avg(R, F, omega_nm, req_nm);
        for (tau = 0 ; tau < NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_nm, req_nm);
            // linear factor
            Linear_num_LSC(omega_nm, R0, R, V0, gamma_array, linear_re, linear_im);
            linear_accum_re[tau] = linear_re;
            linear_accum_im[tau] = linear_im;
            
            MOVEA(R, V, F);
            force_avg(R, F, omega_nm, req_nm);
            MOVEB(V, F);
        }
        
        for (i = 0; i < LEN; i++) {
            if (t_array[i] >= 0) {
                MDlen = static_cast<int>(abs(t_array[i] / DT));
                integral_du[i] = Integrate(du_accum, MDlen, DT); //notice sign of DT
                //new for EFGRL
                mc_re[i] += cos(integral_du[i]/hbar) * linear_accum_re[i] - sin(integral_du[i]/hbar) * linear_accum_im[i];
                mc_im[i] += sin(integral_du[i]/hbar) * linear_accum_re[i] + cos(integral_du[i]/hbar) * linear_accum_im[i];
            }
        }
        
        // ---<<< Backward MD propagation
        for (w = 0; w < n_omega; w++) {
            R[w] = R0[w]; //restore initial condition
            V[w] = V0[w];
        }
        if (t_array[0] < 0) DT = - ABSDT; //for MD propagation direction
        else DT = ABSDT;
        DT2= 0.5 * DT;
        //NMD = NMD_backward;
        
        //dynamics on average surface + wigner sampling
        force_avg(R, F, omega_nm, req_nm);
        for (tau = 0 ; tau< NMD; tau++) {
            //record DU every DT: DU is sum of all frequencies
            du_accum[tau] = DU(R, omega_nm, req_nm);
            //linear factor
            Linear_num_LSC(omega_nm, R0, R, V0, gamma_array, linear_re, linear_im);
            linear_accum_re[tau] = linear_re;
            linear_accum_im[tau] = linear_im;

            MOVEA(R, V, F);
            force_avg(R, F, omega_nm, req_nm);
            MOVEB(V, F);
        }
        
        for (i = 0; i < LEN; i++) {
            if (t_array[i] < 0) {
                MDlen = static_cast<int>(abs(t_array[i]/ DT)); // MDlen should >= 0
                integral_du[i] = Integrate(du_accum, MDlen, DT); //notice sign of DT
                mc_re[i] += cos(integral_du[i]/hbar) * linear_accum_re[i] - sin(integral_du[i]/hbar) * linear_accum_im[i];
                mc_im[i] += sin(integral_du[i]/hbar) * linear_accum_re[i] + cos(integral_du[i]/hbar) * linear_accum_im[i];
            }
        }
        Job_finished(jobdone, j, MCN, startTime);
    }
    
    
    
    //Analyze dynamics on average surface + wigner sampling (LSC)
    for (i = 0; i < LEN; i++) { //Monte Carlo averaging
        mc_re[i] /= MCN;
        mc_im[i] /= MCN;
        corr1[i] = mc_re[i]; //  k(t) re
        corr2[i] = mc_im[i]; //  k(t) im
    }
    
    FFT(-1, mm, corr1, corr2);//notice its inverse FT
    
    for(i = 0; i < nn; i++) { //shift time origin
        corr1_orig[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(-2*pi*i*shift/N);
        corr2_orig[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(-2*pi*i*shift/N);
    }
    
    
    outfile.open((emptystr+"num_LSC_EFGRL_"+idstr+".dat").c_str());
    for (i=0; i<nn/2; i++) outfile << corr1_orig[i]*LEN*DeltaT*DAcoupling*DAcoupling << endl;
    outfile.close();
    outfile.clear();
    

    cout << "-------------------" << endl;
    cout << "DeltaT = " << DeltaT << endl;
    cout << "N = " << LEN << endl;
    cout << "df = " << 1.0/LEN/DeltaT << endl;
    cout << "f_max = " << 0.5/DeltaT << endl;
    cout << "beta = " << beta << endl;
    cout << " eta = " << eta << endl;
    
    cout << "Done: job # " << idstr << endl;



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


//FGR exponent integrand and linear term

void Integrand_exact(double omega, double t, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = (1-cos(omega*t))*Coth;
    im = sin(omega*t);
    return;
}

void Linear_exact(double omega, double t, double req, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = 0.5*hbar/omega*Coth*cos(omega*t) + 0.25*req*req*((1-cos(omega*t))*(1-cos(omega*t)) - sin(omega*t)*sin(omega*t)*Coth*Coth);
    im = -0.5*hbar/omega*sin(omega*t) + 0.5*req*req*Coth*(1-cos(omega*t))*sin(omega*t);
    return;
}

void Integrand_LSC(double omega, double t, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = (1-cos(omega*t))*Coth;
    im = sin(omega*t);
    return;
}

void Linear_LSC(double omega, double t, double req, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = 0.5*hbar/omega* Coth *cos(omega*t) - 0.25*req*req* sin(omega*t)*sin(omega*t)*Coth*Coth;
    im = -0.5*hbar/omega*sin(omega*t) + 0.25*req*req*Coth*(1-cos(omega*t))*sin(omega*t);
    return;
}

void Integrand_CAV(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = sin(omega*t);
    return;
}

void Linear_CAV(double omega, double t, double req, double &re, double &im) {
    re = (beta*hbar*hbar*cos(omega*t) - req*req*sin(omega*t)*sin(omega*t)) / (beta*beta * hbar*hbar *omega*omega);
    im = 0;
    return;
}

void Integrand_CD(double omega, double t, double &re, double &im) {
    re = (1-cos(omega*t))*2/(beta*hbar*omega);
    im = omega*t;
    return;
}

void Linear_CD(double omega, double t, double req, double &re, double &im) {
    re = (beta*hbar*hbar*cos(omega*t) - req*req*sin(omega*t)*sin(omega*t)) / (beta*beta * hbar*hbar *omega*omega);
    im = 0;
    return;
}

void Integrand_W0(double omega, double t, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = omega*omega*t*t*0.5 * Coth;
    im = omega*t;
    return;
}

void Linear_W0(double omega, double t, double req, double &re, double &im) {
    double Coth = 1.0/tanh(beta*hbar*omega*0.5);
    re = ( 0.5*hbar/omega - 0.25*req*req*omega*omega*t*t*Coth)*Coth;
    im = 0;
    return;
}


void Integrand_Marcus(double omega, double t, double &re, double &im) {
    re = omega*t*t/(beta*hbar);
    im = omega*t;
    return;
}
                                                        
void Linear_Marcus(double omega, double t, double req, double &re, double &im) {
  re = 1.0/(beta*omega*omega) - req*req*t*t/(beta*beta*hbar*hbar);
  im = 0;
  return;
}

                                                          
double Integrate(double *data, int n, double dx){
    double I = 0;
    for (int i = 0; i < n; i++) {
        I += data[i];
    }
    I *= dx;
    return I;
}

double Sum(double *data, int n){
    double I = 0;
    for (int i = 0; i < n; i++) {
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
    for (int i = 0; i < N; i++) {
        xx = R[i] + DT*V[i] + DTSQ2*F[i];
        //pbc(xx, yy, zz);
        R[i] = xx;
        V[i] += DT2*F[i];
    }
    return;
}


void MOVEB(double V[], double F[]) {
    //always call MOVEB after call force** to update force F[]
    for (int i = 0; i < N; i++) {
        V[i] += DT2*F[i];
    }
    return;
}

void force_avg(double R[], double F[], double omega[], double req[]) {
    //avg harmonic oscillator potential
    for (int i = 0; i < N; i++) {
        F[i] = - omega[i] * omega[i] * (R[i]- req[i] * 0.5);
    }
    return;
}

void force_donor(double R[], double F[], double omega[], double req[]) {
    //donor harmonic oscillator potential
    for (int i = 0; i < N; i++) {
        F[i] = - omega[i] * omega[i] * R[i];
    }
    return;
}

double DU(double R[], double omega[], double req[]) {
    double du = 0;
    for (int i = 0; i < N; i++) du += req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

double DUi(double R[], double omega[], double req[], int i) {
    double du = 0;
    du = req[i]*omega[i]*omega[i]*R[i]-0.5*req[i]*req[i]*omega[i]*omega[i];
    return du;
}

void Linear_num_LSC(double *omega, double *r0, double *rt, double *p0, double *gm, double &re, double &im) {
    int i;
    re = 0;
    im = 0;
    for (i = 0; i< N ; i++) {
        re += gm[i] *gm[i] * rt[i] * r0[i];
        im -= gm[i] *gm[i] * rt[i] * p0[i] / omega[i] * tanh(beta*hbar*omega[i]*0.5);
    }
    return;
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
