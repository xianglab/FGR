/*test FFT */
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int LEN = 512; //number of t choices
const double DT=0.02;
const double pi=3.14159265358979;

void FFT(int dir, int m, double *x, double *y); //Fast Fourier Transform, 2^m data

int main (int argc, char *argv[]) {
    
    int mm(0), nn(1); // nn = 2^mm is number of (complex) data to FFT
	
	while (nn < LEN ) {
		mm++;
		nn *= 2;
	} //nn is the first 2^m that larger LEN
    
	double *corr1 = new double [nn];//Re
	double *corr2 = new double [nn];//Im
    double d=5.12;
    int i;
    double t;
    double a = 32;
    
    for (i=0; i<LEN; i++) {
        t = i* DT;
        corr1[i] = exp(-a*(t-d)*(t-d)); //Re[f]
        corr2[i] = 0; //Im[f]
    }
    
    FFT(1, mm, corr1, corr2);
    
    ofstream outfile;
    
    outfile.open("FFT-RE.dat");
    cout << "Re[C(t-t_0)] " << endl;
    for (i=LEN/2; i< LEN; i++) {
        outfile << corr1[i] << endl;
    }
    for (i=0; i< LEN/2; i++) {
        outfile << corr1[i] << endl;
    }
    outfile.close();
    
    outfile.open("FFT-IM.dat");
    //cout << " ************* " << endl;
    cout << "Im[C(t-t0)] " << endl;
    for (i=LEN/2; i< LEN; i++) {
        outfile << corr2[i] << endl;
    }
    for (i=0; i< LEN/2; i++) {
        outfile << corr2[i] << endl;
    }
    outfile.close();
    
    double shift = d / DT;
    double N = nn;
    
    double *re_f = new double [nn];//Re
	double *im_f = new double [nn];//Im
    
    for (i=0; i< nn; i++) {
        re_f[i] = corr1[i] * cos(2*pi*i*shift/N) - corr2[i] * sin(2*pi*i*shift/N);
        im_f[i] = corr2[i] * cos(2*pi*i*shift/N) + corr1[i] * sin(2*pi*i*shift/N);
    }
    
    outfile.open("FFT-RE-orig.dat");
    cout << "Re[C(t)] " << endl;
    for (i=LEN/2; i< LEN; i++) {
        outfile << re_f[i] << endl;
    }
    for (i=0; i< LEN/2; i++) {
        outfile << re_f[i] << endl;
    }
    outfile.close();
    
    outfile.open("FFT-IM-orig.dat");
    //cout << " ************* " << endl;
    cout << "Im[C(t)] " << endl;
    for (i=LEN/2; i< LEN; i++) {
        outfile << im_f[i] << endl;
    }
    for (i=0; i< LEN/2; i++) {
        outfile << im_f[i] << endl;
    }
    outfile.close();

    return 0;
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
         \           - i 2 pi k n / N
  X(n) =  >   x(k) e                       = forward transform
         /                                    n=0..N-1
         ---
         k=0
  
  Formula: reverse
              N-1
              ---
           1  \           i 2 pi k n / N
  X(n) =  ---  >   x(k) e                  = forward transform
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
    
	// times STEPS*DT forward FFT (time --> freq)
	/*if (dir == 1) {
		for (i=0; i<n; i++) {
			x[i] *= 1;//DT;
			y[i] *= 1;//DT;
		}
	}*/
    
    // Scaling for inverse transform
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    return;
}
