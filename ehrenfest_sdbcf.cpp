//  ehrenfest_sdbcf.cpp
//
//  " The original Ehrenfest dynamics of 2-level system "
//
//  Compile: g++ -o ehrenfest_sdbcf  ehrenfest_sdbcf.cpp  -std=c++11
//
//  To generate SDBCF < Gamma(a,b,c,d;t) Lambda >_eq
//  defined in Eq. (20) in Shi and Geva, JCP 2004, 120, 10647.
//  call with input a-d values using format: (donor = 0; accep = 1)
//
//  ./ehrenfest_sdbcf [id] [a b c d]
//
//  Created by Xiang Sun on 6/5/16.
//



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>
#include <complex>
using namespace std;



typedef std::complex<double>  Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<double> > Real_Matrix;



//CONSTANTS

const double omega_DA = 2; //0 or 2
const double beta = 1; // temperature
const double Gamma_DA = 0.1; // diabatic coupling (Condon)
const int Ntraj = 20000; // number of independent runs
const int LEN_TRAJ = 2000; // steps in one traj
const double s = 0; // noneq shift
const double DT = 0.01; // time step
const int EPN = 1;
const double dt = DT / EPN;
const int DOFn = 200;
const int DOFe = 2;
const double pi = std::acos(-1.0);
const double hbar = 1;
const double DT2 = 0.5 * DT;
const double DTSQ2 = 0.5 * DT * DT;
const Complex I(0,1);
const int Nobs = 2;



//DECLARE SUBROUTINES
int Initialize_ff(vector<double>& Omega, vector<double>& Req, vector<double>& Shift);

int GenInitialCondition(mt19937& gen, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix& sigma);

int Propagate(int itraj, vector<int>& bcf, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix& sigma, vector<double>& obs, vector<Complex>& sdbcf);

void ForceMF(vector<double>& Omega, vector<double>& Req, vector<double>& R, Complex_Matrix& sigma, vector<double>& Fmf);
double Lambda(vector<double>& R, vector<double>& Req, vector<double>& c);
void Tredof(Complex_Matrix &Ut, vector<int>& bcf, Complex& tre);
void MatrixV(vector<double>& Omega, vector<double>& Req, vector<double>& R, Real_Matrix& V);
void MOVEA(vector<double>& R, vector<double>& P, vector<double>& F);
void MOVEB(vector<double>& P, vector<double>& F);
void MOVEe(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U, Complex_Matrix& sigma);
void MatrixT(Real_Matrix& V, Complex_Matrix& T);
void MatrixU(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U);
void Multiply_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B, Complex_Matrix& C);
void Copy_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B);
void CT_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& Adagger);
int Job_finished(int &jobdone, int count, int total, int startTime);
void OUTPUT_Complex_Matrix(Complex_Matrix& A);
double KE(vector<double>& P);
double Vmf(Complex_Matrix& sigma, Real_Matrix& V);



//MAIN PROGRAM
int main (int argc, char *argv[]) {

    int i,j,k;
    int flag(0);
    int id(0); //job number
    int itraj; //index of trajectory
    int bcf_flag(0);
    vector<int> bcf;//SDBCF indices
    bcf.resize(4,0);
    stringstream ss;
    string emptystr("");
    string idstr("");
    string bcf_ind("");
    
    //condition information
    ss << "b" << beta << "g" << Gamma_DA;
    ss << "_w" << omega_DA << "s" << s;
    idstr += ss.str();
    ss.str("");
    ss.clear();
    
    if (argc > 1) {
        id = atoi(argv[1]);
        ss << "_" << id;
        idstr += ss.str();
        ss.str("");
        ss.clear();
        cout << ">>> Start job id # " << id << endl;
        
        if (argc = 6) {
            bcf[0] = atoi(argv[2]);
            bcf[1] = atoi(argv[3]);
            bcf[2] = atoi(argv[4]);
            bcf[3] = atoi(argv[5]);
            
            for (i=0; i<4; i++)
                if (bcf[i]!=0 && bcf[i]!=1) {
                    cout << "Abort job: Error in SDBCF indices, neither 0 nor 1." << endl;
                    return -1;
                }
            
            bcf_flag = 1;
            
            cout << ">>> Start calculating SDBCF: <Gamma(" << bcf[0] << "," << bcf[1]
            << "," << bcf[2] << "," << bcf[3] <<";t) Lambda>_eq" << endl;
        }
        
    }
    else {
        cout << ">>> ERROR: parameters not enough, abort." << endl;
        return -1;
    }
    
    ss << bcf[0] << bcf[1] << bcf[2] << bcf[3] <<"_";
    bcf_ind = ss.str();
    ss.str("");
    ss.clear();
    
    
    vector<double> R; //nuclear coordinate
    vector<double> P; //nuclear momentum
    vector<double> Omega; //nuclear modes' frequency
    vector<double> Req; //acceptor eq geometry shifts
    vector<double> Shift; //acceptor eq geometry shifts
    vector<double> popd_accum; //the average donor population
    vector<double> E_accum;//the average total energy
    vector<double> bcfre_accum;//sdbcf re
    vector<double> bcfim_accum;//sdbcf im
    vector<double> obs;
    vector<Complex> sdbcf;
    
    Omega.resize(DOFn,0);
    Req.resize(DOFn,0);
    Shift.resize(DOFn,0);
    R.resize(DOFn,0);
    P.resize(DOFn,0);

    popd_accum.resize(LEN_TRAJ,0);
    E_accum.resize(LEN_TRAJ,0);
    bcfre_accum.resize(LEN_TRAJ,0);
    bcfim_accum.resize(LEN_TRAJ,0);
    obs.resize(LEN_TRAJ*Nobs,0);
    sdbcf.resize(LEN_TRAJ, 0.0);
    
    Complex_Matrix sigma(DOFe,vector<Complex>(DOFe,0.0));
    
    
    int jobdone(0);
    int startTime;
    startTime = time(NULL);
    
    random_device rd;
    mt19937 gen; // for random generator
    //gen.seed(id*7+153); //the same seed for debug
    gen.seed(rd());
    
    flag = Initialize_ff(Omega, Req, Shift); //input force field parameters
    
    if (flag != 0) {
        cout << "Terminate due to error in getting force field parameters." << endl;
        return -1;
    }

    
    cout << "    omega_DA = " << omega_DA << endl;
    cout << "    s        = " << s << endl;
    cout << "    beta     = " << beta << endl;
    cout << "    Gamma_DA = " << Gamma_DA << endl;
    cout << "    Ntraj    = " << Ntraj<< endl;
    cout << "    LEN_TRAJ = " << LEN_TRAJ << endl;
    cout << "    DT       = " << DT << endl;

    
    
    //BEGIN proparation of Ntraj Ehrenfest trajectories
    ofstream outfile;
    cout << ">>> Start running Ehrenfest trajectories." << endl;
    
    for (itraj = 0; itraj < Ntraj; itraj++) {
        
        GenInitialCondition(gen, Omega, Req, Shift, R, P, sigma);
     
        //cout << ">> Traj # " << itraj << endl;
        
        Propagate(itraj, bcf, Omega, Req, Shift, R, P, sigma, obs, sdbcf); // run single Ehrenfest trajectory

        for (i = 0 ; i < LEN_TRAJ; i++) {
            
            popd_accum[i] += obs[Nobs*i]; //pop donor
            
            E_accum[i] += obs[Nobs*i + 1]; //total Energy
            
            bcfre_accum[i] += sdbcf[i].real();
            
            bcfim_accum[i] += sdbcf[i].imag();
        }
        
        Job_finished(jobdone, itraj, Ntraj, startTime);
        //cout << "Traj # " << itraj << " done." << endl;
    }
    
    
    outfile.open((emptystr + "SDBCF_re_" + bcf_ind + "Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        bcfre_accum[i] /= Ntraj;
        outfile << bcfre_accum[i] << endl;
    }
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "SDBFC_im_" + bcf_ind + "Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        bcfim_accum[i] /= Ntraj;
        outfile << bcfim_accum[i] << endl;
    }
    outfile.close();
    outfile.clear();
    
    
    cout << ">>> All Traj are done!" << endl;


    return 0;
}





// --------------------------------------------------
// ------------------ SUBROUTINES -------------------
// --------------------------------------------------

int Initialize_ff(vector<double>& Omega, vector<double>& Req, vector<double>& Shift) {
    //for 2-level system
    //read in force field parameters
    
    if (DOFe != 2) cout << "Error: # of electronic states larger than 2" << endl;
    
    int i,j;
    ifstream infile;
    //Omega input
    infile.open("GOA_omega_nm.txt");
    if (!infile.is_open()) {
        cout << "Error: input file Omega cannot open"<< endl;
        return -1;
    }
    if (Omega.size() != DOFn) Omega.resize(DOFn,0);
    for (i = 0; i < DOFn; i++) {
        infile >> Omega[i];
    }
    infile.close();
    infile.clear();
    
    //Req input
    infile.open("GOA_req_nm.txt");
    if (!infile.is_open()) {
        cout << "error: input file Req cannot open"<< endl;
        return -1;
    }
    for (i = 0; i < DOFn; i++) {
        infile >> Req[i];
    }
    infile.close();
    infile.clear();
    
    //Noneq shift input
    infile.open("GOA_shift_nm.txt");
    if (!infile.is_open()) {
        cout << "error: input file Shift cannot open"<< endl;
        return -1;
    }
    for (i = 0; i < DOFn; i++) {
        infile >> Shift[i];
    }
    infile.close();
    infile.clear();
    
    
    cout << ">>> Input force field parameters successfully." << endl;
    
    return 0;
}



int GenInitialCondition(mt19937& gen, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix &sigma) {
    //we label 0 as donor and 1 as acceptor (DOFe=2)

    //[1] initialize nuclear DOF on donor surface (Wigner) at temperature T
    vector<double> sigma_x(DOFn);
    vector<double> sigma_p(DOFn);
    int i,j;
    for (i = 0 ; i < DOFn ; i++) {
        //Wigner distribution:
        sigma_x[i] = sqrt(hbar / (2*Omega[i]*tanh(0.5*beta*hbar*Omega[i])) );
        sigma_p[i] = sigma_x[i] * Omega[i];
    }
    
    normal_distribution<double> normal_dist(0.0,1.0);
    
    for (i = 0 ; i < DOFn ; i++) {
        R[i] = normal_dist(gen) * sigma_x[i] - s * Shift[i];//Wigner initial conf sampling
        P[i] = normal_dist(gen) * sigma_p[i];//Wigner initial momentum sampling
    }
    
    //[2] initialize MM LSC electronic DOF
    for (i = 0 ; i < DOFe ; i++)
        for (j = 0 ; j < DOFe ; j++)
            sigma[i][j] = 0;
    
    sigma[0][0] = 1.0; //initial donor population = 1
    
    return 0;
}



int Propagate(int itraj, vector<int>& bcf, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix& sigma, vector<double>& obs, vector<Complex>& sdbcf){
    //single Ehrenfest trajectory based on initial sampling
    
    int i,j,k;
    double re, im;
    string idstr;
    stringstream ss;
    ss << itraj;
    idstr += ss.str();
    
    //the mean-field force
    vector<double> Fmf(DOFn,0);
    vector<double> c(DOFn,0);
    //potential energy
    Real_Matrix V(DOFe,vector<double>(DOFe,0.0));
    //propagator matrix U
    Complex_Matrix U(DOFe,vector<Complex>(DOFe,0.0));
    //propagator matrix U(t)
    Complex_Matrix Ut(DOFe,vector<Complex>(DOFe,0.0));
    //transformation matrix
    Complex_Matrix T(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp1(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp2(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp3(DOFe,vector<Complex>(DOFe,0.0));
    
    
    Complex Lambda_0w(0.0,0.0); //nuclear part of time 0 SDBCF
    Complex Lambda_t(0.0,0.0); //nuclear part of time t SDBCF
    Complex tre(0.0,0.0); //trace over electronic part of time t SDBCF
    
    for (j = 0; j < DOFn; j++) {//system-bath coupling
        c[j] = 0.5 * Omega[j] * Omega[j] * Req[j];
    }
    
    re = 0;
    im = 0;
    for (j = 0; j < DOFn; j++) {
        im -=  c[j] * P[j] / Omega[j] * tanh(0.5*beta*hbar*Omega[j]);
    }
    
    re = Lambda(R, Req, c);
    Lambda_0w.real(re);
    Lambda_0w.imag(im);

    for (k = 0; k < DOFe; k++) {
        Ut[k][k] = 1;
    }
 
    ForceMF(Omega, Req, R, sigma, Fmf);
    double Lt(0);
    Complex LL(0,0);
    Complex LLe(0,0);
    
    for (i = 0; i < LEN_TRAJ; i++) {
        
        obs[Nobs*i] = sigma[0][0].real();
        //prepare to propagate sigma(t)
        MatrixV(Omega, Req, R, V); // Potential energy
        MatrixT(V, T);
        MatrixU(V, T, U);
        
        Lt = Lambda(R, Req, c);
        Lambda_t.real(Lt);
        LL = Lambda_0w * Lambda_t;
        Tredof(Ut, bcf, tre);
        LLe = LL * tre;
        sdbcf[i] = LLe;
        
        obs[Nobs*i + 1] = KE(P) + Vmf(sigma, V);
        
        //propagate sigma(t) and R(t), P(t)
        MOVEe(V, T, U, sigma);
        
        MOVEA(R, P, Fmf);
        
        ForceMF(Omega, Req, R, sigma, Fmf);
        
        MOVEB(P, Fmf);
        
        Multiply_Complex_Matrix(U, Ut, tmp1);
        
        Copy_Complex_Matrix(tmp1, Ut);
  
    }
    
    return 0;
}



void ForceMF(vector<double>& Omega, vector<double>& Req, vector<double>& R, Complex_Matrix& sigma, vector<double>& Fmf) {
    double fdi, fai;
    for (int i = 0 ; i < DOFn ; i++) {
        fdi = - Omega[i] * Omega[i] * R[i];
        fai = - Omega[i] * Omega[i] * (R[i] - Req[i]);
        Fmf[i] = fdi * sigma[0][0].real() + fai * sigma[1][1].real();
        //Condon approx.
    }
    return;
}



double Lambda(vector<double>& R, vector<double>& Req, vector<double>& c) {
    int i;
    double L(0);
    for (i=0; i< DOFn; i++) {
        L += c[i] * (R[i] - Req[i]*0.5);
    }
    return L;
}



void Tredof(Complex_Matrix &Ut, vector<int>& bcf, Complex& tre) {
    int a,b,c,d;
    a = bcf[0];
    b = bcf[1];
    c = bcf[2];
    d = bcf[3];
    tre = conj(Ut[b][a]) * Ut[c][d];
    return;
}

void MOVEA(vector<double>& R, vector<double>& P, vector<double>& F) {
    for (int i = 0 ; i < DOFn; i++) {
        R[i] += DT*P[i] + DTSQ2*F[i];
        P[i] += DT2*F[i];
    }
    return;
}



void MOVEB(vector<double>& P, vector<double>& F) {
    for (int i = 0 ; i < DOFn ; i++) {
        P[i] += DT2*F[i];
    }
    return;
}



void MOVEe(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U, Complex_Matrix& sigma) {

    Complex_Matrix tmp(DOFe,vector<Complex>(DOFe,0.0));
    Complex_Matrix U_dag(DOFe,vector<Complex>(DOFe,0.0));
    
    Multiply_Complex_Matrix(U, sigma, tmp);
    
    CT_Complex_Matrix(U, U_dag);
    
    Multiply_Complex_Matrix(tmp, U_dag, sigma);
    
    return;
}



void MatrixV(vector<double>& Omega, vector<double>& Req, vector<double>& R, Real_Matrix& V) {
    //two level system
    V[1][0] = V[0][1] = Gamma_DA;
    V[0][0] = omega_DA;
    V[1][1] = 0;
    
    for (int i = 0 ; i < DOFn; i++) {
        V[0][0] += 0.5 * Omega[i] * Omega[i] * R[i] * R[i];
        V[1][1] += 0.5 * Omega[i] * Omega[i] * (R[i] - Req[i]) * (R[i] - Req[i]);
    }
    
    return;
}



void MatrixT(Real_Matrix& V, Complex_Matrix& T) {//two level system
    double t2; // theta/2
    t2 = 0.5 * atan(2*V[0][1] / (V[0][0] - V[1][1]));
    if (t2 < 0) t2 += pi * 0.5;
    T[0][0] = cos(t2);
    T[0][1] = -1 * sin(t2);
    T[1][0] = sin(t2);
    T[1][1] = cos(t2);
    return;
}



void MatrixU(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U) {//two level system
    double Ep, Em;
    double rt, Eav;
    
    rt = sqrt((V[0][0] - V[1][1])*(V[0][0] - V[1][1]) + 4 * V[0][1]* V[1][0] );

    Eav = (V[0][0] + V[1][1]) * 0.5;
    Ep = Eav + rt * 0.5;
    Em = Eav - rt * 0.5;

    Complex_Matrix tmp1(DOFe,vector<Complex>(DOFe,0.0));
    Complex_Matrix tmp2(DOFe,vector<Complex>(DOFe,0.0));
    
    tmp1[0][0] = exp(-I * dt * Ep);
    tmp1[1][1] = exp(-I * dt * Em);
    tmp1[0][1] = tmp1[1][0] = 0;
    
    Multiply_Complex_Matrix(T, tmp1, tmp2);
    CT_Complex_Matrix(T, tmp1);
    Multiply_Complex_Matrix(tmp2, tmp1, U);
    return;
}



void Multiply_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B, Complex_Matrix& C) {
    int i,j,k;
    for (i = 0; i < DOFe; i++)
        for (j = 0; j < DOFe; j++) C[i][j] = 0;
    
    for (i = 0; i < DOFe; i++)
        for (j = 0; j < DOFe; j++) {
            for (k = 0; k < DOFe; k++) C[i][j] += A[i][k] * B[k][j];
        }
    return;
}


void Copy_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B) {
    int i,j;
    for (i = 0; i < DOFe; i++)
        for (j = 0; j < DOFe; j++) B[i][j] = A[i][j];
    
    return;
}


void CT_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& Adagger) {
    int i,j;
    for (i=0; i < DOFe; i++)
        for (j=0; j < DOFe; j++) Adagger[i][j] = conj(A[j][i]);
    return;
}



int Job_finished(int &jobdone, int count, int total, int startTime) {
    int tenpercent;
    int currentTime;
    tenpercent = static_cast<int> (10 * static_cast<double> (count)/ static_cast<double> (total) );
    if ( tenpercent > jobdone ) {
        jobdone = tenpercent;
        currentTime = time(NULL);
        cout << "Finished "<< jobdone <<"0 % : " << count << " out of " << total
        << " trajs.\t Time elapsed " << currentTime - startTime << " sec." << endl;
    }
    return tenpercent;
}



void OUTPUT_Complex_Matrix(Complex_Matrix& A) {
    int i,j;
    for (i=0; i < DOFe; i++) {
        for (j=0; j < DOFe; j++) cout << A[i][j] << "    ";
        cout << endl;
    }
    cout << endl;
    return;
}



double KE(vector<double>& P) {
    double K(0);
    for (int i=0; i < DOFn; i++) {
        K += P[i]*P[i];
    }
    K *= 0.5;
    return K;
}



double Vmf(Complex_Matrix& sigma, Real_Matrix& V) {
    Complex PE(0,0);
    int i,j;
    for (i=0; i < DOFe; i++) {
        for (j=0; j < DOFe; j++) {
            PE += sigma[i][j] * V[j][i];
        }
    }
    
    if (abs(PE.imag()) > 0.00001 ) cout << "Vmf not real" << endl;
    return PE.real();
}






