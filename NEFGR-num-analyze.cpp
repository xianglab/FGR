/* This program generates Non-equilibrium Fermi Golden Rule rate
   in Condon case for Brownian oscillator model with linearized 
   semiclassical method.
   Following multiple runs by ./NEFGR-num #jobid
   To compile: g++ -o NEFGR-num-analyze NEFGR-num-analyze.cpp
   To run: ./NEFGR-num-analyze
   (c) Xiang Sun 2015  */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;

//*********** change parameter *********
string filename = "num_LSC_NEFGR_C_";
const int NFILE = 100;
int MINFILE = 1;
const int MCN = 10000;//50000; //Monte Carlo sample rate
double omega_DA_fix = 0; //fixed omega_DA, with scan tp
double s = 1;    //Noneq. initial shift of parimary mode
const double beta = 5; //3;
const double eta = 0.5;  //3;
double Omega = 0.5;
const int N_BUNCH = 10; //number of bunches for block analysis. (default= 1)
const int BUNCH_FILE = 10;//in each block, # of files. (default = NFILE)
//*********** **************** *********

//double tp_fix = 5; //fixed t' for noneq FGR rate k(t',omega_DA) with scan omega_DA
const double DAcoupling = 0.1;
const double tp_max = 20; //scanning tp option, DeltaTau as step
const double Deltatp = 0.2; //time step for tp (0 ~ tp_max)
const double DeltaTau =0.002; //time step for tau (0~tp griding)

const double pi=3.14159265358979324;
const double hbar = 1;

//*****declare subroutines******
double Integrate(double *data, int n, double dx);
double Sum(double *data, int n);
double** Create_matrix(int row, int col);
void average_error(string filename, int length);


int main (int argc, char *argv[]) {
    
    stringstream ss;
    string emptystr("");
    string idstr("");
    string nameapp("");
    string fnum("");
    string bunchstr("");


    cout << "--------- BEGIN. NEFGR-num-analyze --------" << endl;

    ss.str("");
    nameapp = "";
    ss << "b" << beta << "e" << eta << "_";
    ss << "s" << s << "w" << omega_DA_fix << "_" << MCN << "_";
    nameapp = ss.str();
    
    ss.str("");
    ss << NFILE;
    fnum = ss.str();
    
    int i, j;
    int bunch_index;
    int tp_index_max = static_cast<int> (tp_max / Deltatp);
    int tp_index;
    int M; //time slice for tau = 0 --> tp
    int m; //index of tau
    double tp;
    double tau;
    
    //allocate uncontineous 2D array with different column length
    double **C_re_accum = new double *[tp_index_max];
    double **C_im_accum = new double *[tp_index_max];
    double **C_re = new double *[tp_index_max];
    double **C_im = new double *[tp_index_max];
    for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
        tp = tp_index * Deltatp;
        M = static_cast<int> (tp/DeltaTau);
        C_re_accum[tp_index] = new double [M];
        C_im_accum[tp_index] = new double [M];
        C_re[tp_index] = new double [M];
        C_im[tp_index] = new double [M];
        for (m = 0; m < M; m++) {
            C_re_accum[tp_index][m] = 0; //initialize C_re_accum[][]
            C_im_accum[tp_index][m] = 0; //initialize C_im_accum[][]
        }
    }
    
    double *ktp = new double [tp_index_max]; //transfer rate for BUNCH AVG
    double *ktp_accum = new double [tp_index_max]; //transfer rate accum for AVG
    double *ktp_sq_accum = new double [tp_index_max]; //transfer rate square accum
    for (i=0; i< tp_index_max; i++) {
        ktp_accum[i] = 0;
        ktp_sq_accum[i]=0;
    }
    double *pt = new double [tp_index_max]; // donor state population AVG
    double temp_re;
    double temp_im;

    
    ifstream infile;  //input from file
    ofstream outfile; //output to file  bunch k(tp)
    ofstream outfile1; //output to file  total P(t)
    ofstream outfile2; //output to file  total k(tp)
    ofstream outfile3; //output to file  total k(tp) error
    
    outfile1.open((emptystr+"num_LSC_NEFGR_P_"+nameapp+"AVG"+fnum+".dat").c_str());
    outfile2.open((emptystr+"num_LSC_NEFGR_k_"+nameapp+"AVG"+fnum+".dat").c_str());
    outfile3.open((emptystr+"num_LSC_NEFGR_k_"+nameapp+"ERR"+fnum+".dat").c_str());
    
    
    //BEGIN averaging and error analysis bunch by bunch
    for (bunch_index = 1; bunch_index <= N_BUNCH; bunch_index++)  {
        
        //initialize for accum arrays
        for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
            tp = tp_index * Deltatp;
            M = static_cast<int> (tp/DeltaTau);
            for (m = 0; m < M; m++) {
                C_re_accum[tp_index][m] = 0; //initialize C_re_accum[][]
                C_im_accum[tp_index][m] = 0; //initialize C_im_accum[][]
            }
        }
        
        //import data from files (BUNCH_FILE)
        for (j = 0; j < BUNCH_FILE; j++) {
            ss.str("");
            ss << bunch_index * BUNCH_FILE + MINFILE + j;
            idstr = ss.str();
            
            infile.open((filename+nameapp+idstr+".dat").c_str());
            if (!infile.is_open()) cout << "Error: cannot open data file # " << j << endl;

            for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
                tp = tp_index * Deltatp;
                M = static_cast<int> (tp/DeltaTau);
                for (m = 0; m < M; m++) {
                    infile >> C_re[tp_index][m];
                    infile >> C_im[tp_index][m];
                    C_re_accum[tp_index][m] += C_re[tp_index][m] / BUNCH_FILE;
                    C_im_accum[tp_index][m] += C_im[tp_index][m] / BUNCH_FILE;
                }
            }
            infile.close();
            infile.clear();
        }
        
        //calculate transfer rate k(tp) from C(tp, tau) data (bunch averaging)
        ss.str("");
        ss << bunch_index;
        bunchstr = ss.str();
        outfile.open((emptystr+"num_LSC_NEFGR_k_"+nameapp+"BUNCH"+bunchstr+".dat").c_str());
        for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
            tp = tp_index * Deltatp;
            M = static_cast<int> (tp/DeltaTau);
            ktp[tp_index] = 0;
            //temp_re = temp_im = 0;
            
            for (m = 0; m < M; m++) {
                tau = m * DeltaTau;
                temp_re = C_re_accum[tp_index][m];
                temp_im = C_im_accum[tp_index][m];
                ktp[tp_index] += temp_re * cos(omega_DA_fix*tau) - temp_im * sin(omega_DA_fix*tau);
            }
            ktp[tp_index] *= DeltaTau * 2.0 * DAcoupling*DAcoupling;
            outfile << ktp[tp_index] << endl;
            //accumulate for total average and error
            ktp_accum[tp_index] += ktp[tp_index];
            ktp_sq_accum[tp_index] += ktp[tp_index] * ktp[tp_index];
            //population += ktp[tp_index] * Deltatp; // population of donor state
            //outfile1 << exp(-1.0 * population) << endl; //P(t) = exp(- int_0^t dt' k(t'))
        }
        outfile.close();
        outfile.clear();
        
    }
    double population(0); //population of Donor state
    double error;
    //total averaging and output
    for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
        tp = tp_index * Deltatp;
        
        ktp_accum[tp_index] /= N_BUNCH;
        ktp_sq_accum[tp_index] /= N_BUNCH;
        
        population += ktp_accum[tp_index] * Deltatp; // population of donor state
        outfile1 << exp(-1.0 * population) << endl; //P(t) = exp(- int_0^t dt' k(t'))
        
        outfile2 << ktp_accum[tp_index] << endl;//total avg k(tp)
        
        error = (ktp_sq_accum[tp_index] - ktp_accum[tp_index] * ktp_accum[tp_index]) / sqrt(N_BUNCH);
        
        outfile3 << error << endl; //error bar of k(tp)

    }
    
    
    /*
     //import data from files
     for (j = 0; j < NFILE; j++) {
     ss.str("");
     ss << j + MINFILE;
     idstr = ss.str();
     
     infile.open((filename+nameapp+idstr+".dat").c_str());
     if (!infile.is_open()) cout << "Error: cannot open data file # " << j << endl;
     
     for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
     tp = tp_index * Deltatp;
     M = static_cast<int> (tp/DeltaTau);
     for (m = 0; m < M; m++) {
     infile >> C_re[tp_index][m];
     infile >> C_im[tp_index][m];
     C_re_accum[tp_index][m] += C_re[tp_index][m] / NFILE;
     C_im_accum[tp_index][m] += C_im[tp_index][m] / NFILE;
     }
     }
     infile.close();
     infile.clear();
     }
     
     
     //calculate transfer rate from C(tp, tau) data
     outfile1.open((emptystr+"num_LSC_NEFGR_P_"+nameapp+"AVG"+fnum+".dat").c_str());
     outfile.open((emptystr+"num_LSC_NEFGR_k_"+nameapp+"AVG"+fnum+".dat").c_str());
     
     double *ktp = new double [tp_index_max]; //transfer rate
     double *pt = new double [tp_index_max]; // donor state population
     double temp_re;
     double temp_im;
     double population(0); //population of Donor state
     
     for (tp_index = 0; tp_index < tp_index_max; tp_index++) {
     tp = tp_index * Deltatp;
     M = static_cast<int> (tp/DeltaTau);
     ktp[tp_index] = 0;
     //temp_re = temp_im = 0;
     
     for (m = 0; m < M; m++) {
     tau = m * DeltaTau;
     temp_re = C_re_accum[tp_index][m];
     temp_im = C_im_accum[tp_index][m];
     ktp[tp_index] += temp_re * cos(omega_DA_fix*tau) - temp_im * sin(omega_DA_fix*tau);
     }
     ktp[tp_index] *= DeltaTau * 2.0 * DAcoupling*DAcoupling;
     outfile << ktp[tp_index] << endl;
     population += ktp[tp_index] * Deltatp; // population of donor state
     //outfile1 << 1 - population << endl; //P(t) = 1 - int_0^t dt' k(t')
     outfile1 << exp(-1.0 * population) << endl; //P(t) = exp(- int_0^t dt' k(t'))
     }
     
     outfile.close();
     outfile.clear();
     outfile1.close();
     outfile1.clear();
     */
    
    cout << "   Parameters: " << endl;
    cout << "       beta = " << beta << endl;
    cout << "        eta = " << eta << endl;
    cout << "       fix omega_DA = " << omega_DA_fix << endl;
    cout << "       initial shift s = " << s << endl;
    cout << "   N_BUNCH = " << N_BUNCH << endl;
    cout << "   Averaged " << NFILE << " files." << endl;
    cout << "--------- DONE. NEFGR-num-analyze --------" << endl;
    return 0;
}



/********* SUBROUTINES *************/

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



double** Create_matrix(int row, int col) {//allocate continuous memory for 2D matrix
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (int i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    return matrix;
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
        ss << j + MINFILE;
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
















