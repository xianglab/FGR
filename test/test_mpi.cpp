#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "mpi.h"
using namespace std;

int main (int argc, char *argv[]) {
	
	//******************************************************* 
	//*         Initiate MPI parallel enviornment			*
	//******************************************************* 	
	int err(0);
	int myrank(0), numprocs(0);
	int dest(1);
	int namelen,version, subversion;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double start_time(0), end_time(0);
	MPI_Request request;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);
	cout << "Node: " << processor_name <<"  Processor #" << myrank <<" of "<<numprocs<<endl;
	
	//============ROOT processor (rank=0)========================================================
    if (myrank == 0) {
        cout << "ROOT processor is here" << endl;
    }

    //===============SLAVE processors  (rank=1,2,3,...,numprocs-1)	=====================================
    else {
    }

        MPI_Finalize();
        return (0);
    }   
