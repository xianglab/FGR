#include <iostream>
#include <cstdlib> //atoi:cstring to int
#include <sstream>

using namespace std;

int main (int argc, char *argv[]) {
    int i;
    int j,k(0);
    //istringstream iss;
    stringstream ss;
    cout << "# of argument: " << argc-1 << endl;
    if (argc > 1) {
        for (i=1; i<argc; i++) {
            //iss.str(argv[i]) ;
            //iss >> j;
            //if (iss.fail()) cout << "error " << endl;
            //iss.clear();
            
            ss << argv[i] ;
            ss >> j;
            ss.clear();
            //j = atoi(argv[i]);
            k += j;
        }
        cout << "sum = " << k << endl;
    }
    return 0;
}
