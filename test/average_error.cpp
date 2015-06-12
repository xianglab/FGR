/* This code calculates average and error bars of a serious input files
 (c) Xiang Sun 2014
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;

const int NFILE = 100;
const int MIN = 1;

const double omega_max = 2.5; //20; //20 for ohmic, 2.5 for gaussian

const double d_omega = 0.025;
const int n_omega = static_cast<int>(omega_max / d_omega);
const int N = n_omega; //number of degrees of freedom
const int LEN = 1024; //number of time steps in FFT
const int STEPS = 50; //for correlation function skip steps
const double DeltaT=0.3; //FFT time sampling interval
const double T0= -153.6; //FFT time sampling starting point
double DT=0.01; //MD time step

void average_error(string filename, int length);

int main (int argc, char *argv[]) {
    
    int id; //jobid
    int i,j,k;
    stringstream ss;
    string filename;
    string nameapp = "Gau_b4_5e4_";
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
    
    int NMD;
    int NMD_forward = static_cast<int>(abs(t_array[LEN-1]/DT));
    int NMD_backward = static_cast<int>(abs(t_array[0]/DT));
    if(NMD_forward > NMD_backward) NMD = NMD_forward; //get max of NMD_forward and NMD_backward
    else NMD = NMD_backward;
    
    //time domain real part
    
    filename = "2cumu_ana_t_re_"+nameapp;
    average_error(filename,LEN);
    
    filename = "2cui_t_re_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_re_"+nameapp;
    average_error(filename,LEN);

    filename = "nlsc_t_re_D_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_re_CL_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_re_CL_D_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_re_inh_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_re_inh_CL_"+nameapp;
    average_error(filename,LEN);
    
    filename = "QM_t_re_"+nameapp;
    average_error(filename,LEN);

    //time domain imaginary part
    filename = "2cumu_ana_t_im_"+nameapp;
    average_error(filename,LEN);
    
    filename = "2cui_t_im_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_D_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_CL_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_CL_D_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_inh_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_t_im_inh_CL_"+nameapp;
    average_error(filename,LEN);
    
    filename = "QM_t_im_"+nameapp;
    average_error(filename,LEN);

    //frequency domain real part
    filename = "2cumu_ana_re_"+nameapp;
    average_error(filename,LEN);
    
    filename = "2cui_re_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_re_"+nameapp;
    average_error(filename,NN);

    filename = "nlsc_re_D_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_re_CL_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_re_CL_D_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_re_inh_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_re_inh_CL_"+nameapp;
    average_error(filename,NN);
    
    filename = "QM_re_"+nameapp;
    average_error(filename,NN);

    //frequency domain imaginary part
    filename = "2cumu_ana_im_"+nameapp;
    average_error(filename,LEN);
    
    filename = "2cui_im_"+nameapp;
    average_error(filename,LEN);
    
    filename = "nlsc_im_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_im_D_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_im_CL_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_im_CL_D_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_im_inh_"+nameapp;
    average_error(filename,NN);
    
    filename = "nlsc_im_inh_CL_"+nameapp;
    average_error(filename,NN);
    
    filename = "QM_im_"+nameapp;
    average_error(filename,NN);

    //Cuu time correlation function
    filename = "cuu_"+nameapp;
    average_error(filename, NMD);
    
    filename = "cuu_D_"+nameapp;
    average_error(filename, NMD);
    
    filename = "cuu_eq_"+nameapp;
    average_error(filename, NMD);
    
    filename = "cuu_eq_D_"+nameapp;
    average_error(filename, NMD);
    
    
    cout << "Averaged " << NFILE << " files." << endl;
    
    return 0;
}


// ***************** SUBROUTINE *******************
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



