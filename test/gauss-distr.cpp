#include <cmath>
#include <cstdlib>      /* srand, rand */
#include <ctime>        /* time */
#include <iostream>
#include <fstream>
#include "r_1279.h"  // r1279 random number generator
using namespace std;

double ran2(long *idum);
double GAUSS(long *seed);
double f(double x); //integrad function

int main () {
    int i,j,k;
    double s,r;
    int n;
    n=0;
    long seed;
    s=0;
    
    //seed = - static_cast<long> (time(0));
    seed 	= seedgen();	/* have seedgen compute a random seed */
    setr1279(seed);		/* seed the genertor */

    r1279();
    
    for (i=0; i< 100000000; i++) {
        //s = ran2(&seed);
        //cout << s << endl;
        //cout << GAUSS(&seed) << endl;
        //cout << r1279() << endl;
        s += f(GAUSS(&seed));
        //s += f(r1279());
        n++;
    }
    
    //srand( static_cast<unsigned int> (time(0))); // random seed
    /*
    for (i=0; i< 100000; i++) {
        s = rand()/(double)RAND_MAX;
        cout << s << endl;
    }*/
    cout << " Integral = " << s / n << endl;
    
    return 0;
}

//GENERATE UNIFORM RANDOM NUMBER
double ran2(long *idum) {
    const int IM1=2147483563;
    const int IM2=2147483399;
    const int IA1=40014;
    const int IA2=40692;
    const int IQ1=53668;
    const int IQ2=52774;
    const int IR1=12211;
    const int IR2=3791;
    const int IMM1=(IM1-1);
    const int NTAB=32;
    const int NDIV=(1+IMM1/NTAB);
    const double EPS=1.2e-7;
    const double RNMX=(1.0-EPS);
    const double AM=1.0/(double)IM1;
    
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float tem;
    if (*idum<=0) {
        if (-(*idum)<1) *idum = 1;
		else *idum = -(*idum);
        idum2 = (*idum);
        for (j=NTAB+7; j>0; j--) {
            k = (*idum)/IQ1;
            *idum = IA1*(*idum - k*IQ1) - k*IR1;
            if (*idum<0) *idum += IM1;
            if (j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum = IA1*(*idum-k*IQ1) - k*IR1;
    
    if (*idum<0) *idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2) - k*IR2;
    
    if (idum2<0) idum2 += IM2;
    j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy<1) iy += IMM1;
    if ((tem=AM*iy)>RNMX) return (RNMX);
    else return (tem);
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

double f(double x) {
    return x*x;
}
