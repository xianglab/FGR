//  mf_orig.cpp
//
//  " The original Ehrenfest dynamics of 2-level system "
//
//  Compile: g++ -std=c++11 mf_orig.cpp -o mf_orig
//
//  Created by Xiang Sun on 5/19/16.
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

const double omega_DA = 2;
const int bldim = 3;//3;
const int gldim = 3;
const int sldim = 4;
double beta_list[bldim] = {1, 2, 5}; //{1, 2, 5};
double gamma_list[gldim] = {0.05, 0.1, 0.2};
double s_list[sldim] = {-1.0, 1.0, 3.0, 5.0};
double beta = 1; // 1/temperature
double Gamma_DA = 0.1; //diabatic coupling (Condon)
double s = 1; // noneq shift

const int Ntraj = 10000;//0000;// 10000; // number of independent runs
const int LEN_TRAJ = 1000;//10000;//10000;//10000; //nuclear steps in one traj
const double DT = 0.01;//0.001; // nuclear time step
const int EPN = 1; // number of electronic steps per nuclear step
const double dt = DT / EPN; // electronic time step
const int DOFn = 200;
const int DOFe = 2; //two level system
const double pi = std::acos(-1.0);//3.14159265358979;
const double hbar = 1;
const double DT2 = 0.5 * DT;
const double DTSQ2 = 0.5 * DT * DT;
const Complex I(0,1); //the imaginary I
const int Nobs = 5; //number of data per time step for observable



//DECLARE SUBROUTINES
int Initialize_ff(vector<double>& Omega, vector<double>& Req, vector<double>& Shift);

int GenInitialCondition(mt19937& gen, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix& sigma);

int Propagate(int itraj, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix& sigma, vector<double>& obs);

void ForceMF(vector<double>& Omega, vector<double>& Req, vector<double>& R, Complex_Matrix& sigma, vector<double>& Fmf);
void MatrixV(vector<double>& Omega, vector<double>& Req, vector<double>& R, Real_Matrix& V);
void MOVEA(vector<double>& R, vector<double>& P, vector<double>& F);
void MOVEB(vector<double>& P, vector<double>& F);
void MOVEe(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U, Complex_Matrix& sigma);
void MatrixT(Real_Matrix& V, Complex_Matrix& T);
void MatrixU(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U);
void Multiply_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B, Complex_Matrix& C);
void CT_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& Adagger);
int Job_finished(int &jobdone, int count, int total, int startTime);
void OUTPUT_Complex_Matrix(Complex_Matrix& A);
void OUTPUT_Real_Matrix(Real_Matrix& A);
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
    

    
    if (argc > 1) {
        id = atoi(argv[1]);
        ss << "_" << id;
        idstr += ss.str();
        ss.str("");
        ss.clear();
        cout << ">>> Start job id # " << id << endl;
    }

    vector<double> R; //nuclear coordinate
    vector<double> P; //nuclear momentum
    vector<double> Omega; //nuclear modes' frequency
    vector<double> Req; //acceptor eq geometry shifts
    vector<double> Shift; //acceptor eq geometry shifts
    vector<double> sigma_accum; //the average reduced density matrix
    vector<double> E_accum;//the average total energy
    vector<double> obs;
    
    Omega.resize(DOFn,0);
    Req.resize(DOFn,0);
    Shift.resize(DOFn,0);
    R.resize(DOFn,0);
    P.resize(DOFn,0);

    sigma_accum.resize(LEN_TRAJ * 4 ,0);//pop_donor, sigma_x, sigma_y, sigma_z
    E_accum.resize(LEN_TRAJ,0);
    obs.resize(LEN_TRAJ * Nobs, 0);
    
    //reduced electronic density matrix (2 x 2) complex matrix
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

   
    //BEGIN proparation of Ntraj Ehrenfest trajectories
    ofstream outfile;
    ofstream outfile1;
    cout << ">>> Start running Ehrenfest trajectories." << endl;
    
    //outfile.open((emptystr + "initial_condition_R.dat").c_str());
    //outfile1.open((emptystr + "initial_condition_P.dat").c_str());
    
    
    int beta_index(0);
    int gamma_index(0);
    int s_index(0);

    //BEGIN loop through different conditions
    int case_count(0);
    for (beta_index = 0; beta_index < bldim; beta_index++)
        for (gamma_index = 0; gamma_index < gldim; gamma_index++)
            for (s_index = 0; s_index < sldim; s_index++)
    {
        beta = beta_list[beta_index];
        Gamma_DA = gamma_list[gamma_index];
        s = s_list[s_index];
    
        //condition information
        ss << "b" << beta << "g" << Gamma_DA;
        ss << "_w" << omega_DA << "s" << s;
        //ss << "_" << Ntraj;
        idstr = ss.str();
        ss.str("");
        ss.clear();

        for (i = 0 ; i < LEN_TRAJ; i++) {
            
            sigma_accum[4*i] = 0 ; //pop_donor
            
            sigma_accum[4*i+1] = 0; //sigma_x
            
            sigma_accum[4*i+2] = 0; // sigma_y
            
            sigma_accum[4*i+3] = 0; // sigma_z
            
            E_accum[i] = 0; //total Energy
        }
        
        
        
    for (itraj = 0; itraj < Ntraj; itraj++) {
        
        GenInitialCondition(gen, Omega, Req, Shift, R, P, sigma);
     
        //cout << ">> Traj # " << itraj << endl;
        
        //output initial condition
        //for (i = 0 ; i < DOFn ; i++) {
        //    outfile << R[i] << endl;
        //    outfile1 << P[i] << endl;
        //}
        
        
        Propagate(itraj, Omega, Req, Shift, R, P, sigma, obs); // run single Ehrenfest trajectory

        for (i = 0 ; i < LEN_TRAJ; i++) {
            
            sigma_accum[4*i] += obs[Nobs*i]; //pop_donor
            
            sigma_accum[4*i+1] += obs[Nobs*i+1]; //sigma_x
            
            sigma_accum[4*i+2] += obs[Nobs*i+2]; // sigma_y
            
            sigma_accum[4*i+3] += obs[Nobs*i+3]; // sigma_z
            
            E_accum[i] += obs[Nobs*i + 4]; //total Energy
        }
        
        //Job_finished(jobdone, itraj, Ntraj, startTime);
        
        //cout << "Traj # " << itraj << " done." << endl;
    }
    
    //outfile.close();
    //outfile.clear();
    //outfile1.close();
    //outfile1.clear();
    
    
    outfile.open((emptystr + "PopD_Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        sigma_accum[4*i] /= Ntraj;
        outfile << setprecision(8) << sigma_accum[4*i] << endl;
    }
    outfile.close();
    outfile.clear();
    
    /*
    outfile.open((emptystr + "sigma_x_Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        sigma_accum[4*i+1] /= Ntraj;
        outfile << setprecision(8) << sigma_accum[4*i+1] << endl;
    }
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "sigma_y_Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        sigma_accum[4*i+2] /= Ntraj;
        outfile << setprecision(8) << sigma_accum[4*i+2] << endl;
    }
    outfile.close();
    outfile.clear();
    
    outfile.open((emptystr + "sigma_z_Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        sigma_accum[4*i+3] /= Ntraj;
        outfile << setprecision(8) << sigma_accum[4*i+3] << endl;
    }
    outfile.close();
    outfile.clear();
    
    
    outfile.open((emptystr + "energy_Ehr_" + idstr + ".dat").c_str());
    for (i = 0 ; i < LEN_TRAJ; i++) {
        E_accum[i] /= Ntraj;
        outfile << setprecision(8) << E_accum[i] << endl;
    }
    outfile.close();
    outfile.clear();
    */
    
        
    //-------------- Summary ----------------
    case_count++;
        
    cout << ">>> CASE # " << case_count <<  " done:" << endl;
    cout << "    beta     = " << beta << endl;
    cout << "    Gamma_DA = " << Gamma_DA << endl;
    cout << "    s(shift) = " << s << endl;

    }

    cout << ">>> ---- Universal Parameters ---- " << endl;
    cout << "    omega_DA = " << omega_DA << endl;
    cout << "    Ntraj    = " << Ntraj<< endl;
    cout << "    LEN_TRAJ = " << LEN_TRAJ << endl;
    cout << "    DT       = " << DT << endl;
    cout << ">>> Finished all Ehrenfest trajectories! ----" << endl;
    
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
        //classical distribution:
        /*
         sigma_x[i] = 1/(Omega[i]*sqrt(beta));
         sigma_p[i] = 1/sqrt(beta);
         */
    }
    //normal distribution random number generator c++11
    /*
    random_device rd;
    mt19937 gen(rd());
    */
    //mt19937 gen;
    //gen.seed(seed);
    
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



int Propagate(int itraj, vector<double>& Omega, vector<double>& Req, vector<double>& Shift, vector<double>& R, vector<double>& P, Complex_Matrix &sigma, vector<double>& obs) {
    //single Ehrenfest trajectory based on initial sampling
    
    int i,j,k;
    
    //the mean-field force
    vector<double> Fmf(DOFn,0);
    //potential energy matrix (2 x 2) real matrix
    Real_Matrix V(DOFe,vector<double>(DOFe,0.0));
    //propagator matrix U = exp{-i dt /hbar * V[R(t)]}
    Complex_Matrix U(DOFe,vector<Complex>(DOFe,0.0));
    //transformation matrix
    Complex_Matrix T(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp1(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp2(DOFe,vector<Complex>(DOFe,0.0));
    
    Complex_Matrix tmp3(DOFe,vector<Complex>(DOFe,0.0));
    
    
    
    string idstr;
    stringstream ss;
    ss << itraj;
    idstr += ss.str();
    
    ofstream outfile1;
    //outfile1.open((idstr+"_RPF[50](t).dat").c_str());
    
    //propagate nuclear DOF
    ForceMF(Omega, Req, R, sigma, Fmf); //Fmf(t) using R(t), sigma(t)
    
    /*
    cout << "Force at time 0: " << endl;
    //test for time 0 force vector
    for (j=0; j< DOFn; j++) {
        cout << Fmf[j]<< endl;
    }
    */
    

    
    for (i = 0; i < LEN_TRAJ; i++) {
        
        //record observable of the current configuration R(t) and sigma(t)
        obs[Nobs*i] = sigma[0][0].real(); //sigma[0][0].real(); or norm(sigma[0][0])
        obs[Nobs*i+1] = sigma[0][1].real() + sigma[1][0].real(); //sigma_x
        obs[Nobs*i+2] = sigma[0][1].imag() - sigma[1][0].imag(); //sigma_y
        obs[Nobs*i+3] = sigma[0][0].real() - sigma[1][1].real(); //sigma_z

        //test conservaton of total population
        //cout << "sigma_matrix = " << endl;
        //OUTPUT_Complex_Matrix(sigma);
        //cout << i << " step: total Pop = " << sigma[0][0]+sigma[1][1] << endl;
       
        
        //test for 1 traj: mode 50
        //outfile1 << R[50] << "  "<< P[50] << "  " << Fmf[50] << endl;
        
        //prepare to propagate sigma(t)
        MatrixV(Omega, Req, R, V); // Potential energy matrix V(t) <--- R(t)
        MatrixT(V, T);
        MatrixU(V, T, U);
        
        //if (i<11) {
        //    cout << "Time step " << i << ": V matrix = " << endl;
        //    OUTPUT_Real_Matrix(V);
        //}
        
        
        //test: output total energy, should be conserved
        obs[Nobs*i + 4] = KE(P) + Vmf(sigma, V);
        
        //test: U matrix diagonalize
        /*
        for (j=0;j<2;j++)
            for (k=0; k<2; k++)
                    tmp2[j][k] = V[j][k];
        
        CT_Complex_Matrix(T, tmp1);
        Multiply_Complex_Matrix(tmp1 ,tmp2, tmp3);
        Multiply_Complex_Matrix(tmp3 ,T, tmp1);

        cout << "Tdag * V * T =  diag{E+, E-}:" << endl;
        OUTPUT_Complex_Matrix(tmp1);
        */
        
        //propagate electronic DOF for EPN time steps
        //for (j = 0 ; j < EPN ; j++) {
            //sigma(t+dt) <--- sigma(t), R(t)
            //
        //}
        
        MOVEe(V, T, U, sigma);
        
        MOVEA(R, P, Fmf);//R(t+DT), P(t+DT/2) <--- R(t), P(t), Fmf(t)
        
        ForceMF(Omega, Req, R, sigma, Fmf); //Fmf(t+DT) using R(t+DT), sigma(t+DT)
        
        MOVEB(P, Fmf); // P(t+DT) <--- P(t+DT/2), Fmf(t+DT)
        
  
    }
    
    return 0;
}



void ForceMF(vector<double>& Omega, vector<double>& Req, vector<double>& R, Complex_Matrix& sigma, vector<double>& Fmf) {
    double fdi, fai;
    for (int i = 0 ; i < DOFn ; i++) {
        fdi = (- Omega[i]) * Omega[i] * R[i];
        fai = (- Omega[i]) * Omega[i] * (R[i] - Req[i]);
        Fmf[i] = fdi * sigma[0][0].real() + fai * sigma[1][1].real();
        //Condon approx.
    }
    return;
}



void MOVEA(vector<double>& R, vector<double>& P, vector<double>& F) {
    //VELOCITY VERLET ALGORITHM Part A
    //R(t+DT) = R(t) + V(t)*DT + F(t) * DT^2/2
    //P(t+DT/2) = P(t) + F(t) * DT / 2
    for (int i = 0 ; i < DOFn; i++) {
        R[i] += DT*P[i] + DTSQ2*F[i];
        P[i] += DT2*F[i];
    }
    return;
}



void MOVEB(vector<double>& P, vector<double>& F) {
    //VELOCITY VERLET ALGORITHM Part B
    //always call ForceMF to update force F(t+DT), before call MOVEB
    //P(t+DT) = P(t+DT/2) + F(t+DT) * DT /2
    for (int i = 0 ; i < DOFn ; i++) {
        P[i] += DT2*F[i];
    }
    return;
}



void MOVEe(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U, Complex_Matrix& sigma) {

    Complex_Matrix tmp(DOFe,vector<Complex>(DOFe,0.0));
    Complex_Matrix U_dag(DOFe,vector<Complex>(DOFe,0.0));
    
    Multiply_Complex_Matrix(U, sigma, tmp); //tmp = U * sigma(t)
    
    CT_Complex_Matrix(U, U_dag); // U_dag = U^dagger
    
    Multiply_Complex_Matrix(tmp, U_dag, sigma); // sigma(t+dt) = tmp * U^dagger
    
    return;
}



void MatrixV(vector<double>& Omega, vector<double>& Req, vector<double>& R, Real_Matrix& V) {
    //two level system
    V[1][0] = V[0][1] = Gamma_DA;
    V[0][0] = omega_DA;
    V[1][1] = 0;
    
    //test: one mode only i=11
    //int i=11;
    //V[0][0] += 0.5 * Omega[i] * Omega[i] * R[i] * R[i];
    //V[1][1] += 0.5 * Omega[i] * Omega[i] * (R[i] - Req[i]) * (R[i] - Req[i]);
    
    
    for (int i = 0 ; i < DOFn; i++) {
        V[0][0] += 0.5 * Omega[i] * Omega[i] * R[i] * R[i];
        V[1][1] += 0.5 * Omega[i] * Omega[i] * (R[i] - Req[i]) * (R[i] - Req[i]);
    }
    
    return;
}



void MatrixT(Real_Matrix& V, Complex_Matrix& T) {//two level system
    double t2; // theta/2
    //Complex eiphi; //e^{i phi} assume V[0][1]=V[1][0]=real
    
    t2 = 0.5 * atan(2*V[0][1] / (V[0][0] - V[1][1]));
    
    if (t2 < 0) t2 += pi * 0.5;
    
    //eiphi = 1;//V[0][1] / abs(V[0][1]);
    
    T[0][0] = cos(t2); // * conj(eiphi);
    T[0][1] = -1 * sin(t2); //* conj(eiphi);
    T[1][0] = sin(t2); // * eiphi;
    T[1][1] = cos(t2); // * eiphi;
    
    return;
}



void MatrixU(Real_Matrix& V, Complex_Matrix& T, Complex_Matrix& U) {//two level system
    
    double Ep, Em; //E+ / E-
    double rt, Eav;
    
    rt = sqrt((V[0][0] - V[1][1])*(V[0][0] - V[1][1]) + 4 * V[0][1]* V[1][0] );
    //norm(Complex z)= |z|^2
    Eav = (V[0][0] + V[1][1]) * 0.5;
    Ep = Eav + rt * 0.5;
    Em = Eav - rt * 0.5;
    
    //working matrix
    Complex_Matrix tmp1(DOFe,vector<Complex>(DOFe,0.0));
    Complex_Matrix tmp2(DOFe,vector<Complex>(DOFe,0.0));
    
    tmp1[0][0] = exp(-I * dt * Ep);
    tmp1[1][1] = exp(-I * dt * Em);
    tmp1[0][1] = tmp1[1][0] = 0;
    
    Multiply_Complex_Matrix(T, tmp1, tmp2); // tmp2 = T * diag{e^{-i dt E+}, e^{-i dt E-} }
    
    CT_Complex_Matrix(T, tmp1); // tmp1 = T^dagger
    
    Multiply_Complex_Matrix(tmp2, tmp1, U); // U = tmp2 * T^dagger
    
    return;
}



void Multiply_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& B, Complex_Matrix& C) {
    int i,j,k;
    for (i=0; i< DOFe; i++)
        for (j=0; j< DOFe; j++) C[i][j] = 0;
    
    for (i=0; i< DOFe; i++)
        for (j=0; j< DOFe; j++) {
            for (k=0; k< DOFe; k++) C[i][j] += A[i][k] * B[k][j];
        }
    return;
}



void CT_Complex_Matrix(Complex_Matrix& A, Complex_Matrix& Adagger) {
    //conjugate transpose of complex matrix
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
    //conjugate transpose of complex matrix
    int i,j;
    for (i=0; i < DOFe; i++) {
        for (j=0; j < DOFe; j++) cout << A[i][j] << "    ";
        cout << endl;
    }
    cout << endl;
    return;
}


void OUTPUT_Real_Matrix(Real_Matrix& A) {
    //conjugate transpose of real matrix
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
            //PE += sigma[j][j] * V[j][j]; //without coherence
            PE += sigma[i][j] * V[j][i];
        }
    }
    
    if (abs(PE.imag()) > 0.00001 ) cout << "Vmf not real" << endl;
    return PE.real();
}






