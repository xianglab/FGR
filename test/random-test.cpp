/* ************ COPYRIGHT 2011, H. G. Katzgraber. ********************* */
/*									*/
/* EX 1:        Learn to use the r1279 RNG and seed generator           */
/*									*/
/*		Computes 10^N random numbers and calculates average	*/
/*		and deviation from exact value, which should be 0.5	*/
/*		Feel free to change N and see how the estimate changes	*/
/* INSTR:	Compile with Makefile via 'make', execute with ./runme	*/
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "r_1279.h"
using namespace std;

#define N	6           	/* number of 10^N random nums    	*/


int main(){
    
    int		k,m;
    double 	sum;
    long 	seed;
    
    seed 	= seedgen();	/* have seedgen compute a random seed 	*/
    setr1279(seed);		/* seed the genertor			*/
    
    m 		= pow(10,N);	/* compute number of random numbers 	*/
    
    sum 	= 0.0;
    
    for(k = 1; k <= m; k++){
        sum += r1279();
    }
    
    sum /= ((double) m);
    
    cout << "average = " << sum << ", deviation from exact = " <<sum - 0.5 << endl;
    
    return 0;
}