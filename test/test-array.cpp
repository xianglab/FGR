#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
void func (double *data, int n);

int main (int argc, char *argv[]) {
    int i;
    double x[100];
    for (i=0;i<100;i++) x[i]=0.1*i;
    
    func(&(x[10]), 5);

    return 0;
}

void func (double *data, int n) {
    int j;
    for (j=0; j<n;j++) cout << data[j] << endl;
    return;
}