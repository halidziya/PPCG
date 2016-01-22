#include "Matrix.h"
#include "util.h"
#include "DebugUtils.h"
#include "DataSet.h"
#include "ThreadPool.h"
#include <iostream>
#include <fstream>
#include <list>
#include "IWishart.h"
#include "Algorithms.h"
using namespace std;


int MAX_SWEEP = 500;
int BURNIN = 100;
int SAMPLE = (MAX_SWEEP - BURNIN) / 10; // Default value is 10 sample + 1 post burnin
char* result_dir = "results";

// Variables
int d, n;



int main(int argc, char** argv)
{


	if (argc < 4) {
		printf("Usage igmm datafilename priorfilename configfilename");
		exit(-1);
	}


	char* datafile(argv[1]);
	char* priorfile(argv[2]);
	char* configfile(argv[3]);
	
	

	if (argc>4)
		MAX_SWEEP = atoi(argv[4]);
	if (argc>5)
		BURNIN = atoi(argv[5]);
	if (argc > 6)
		result_dir = argv[6];
	SAMPLE = (MAX_SWEEP - BURNIN) / 10; // Default value
	if (argc>7)
		SAMPLE = atoi(argv[7]);

	step();
	
	string ss(result_dir);
	printf("size of int %d", sizeof(int));
	printf("Reading...\n");
	DataSet ds(datafile, priorfile, configfile);
	ofstream nsampleslog(ss.append("nsamples.log"));
	d = Global::d;
	n = ds.n;
	
	Vector priormean(d);
	Matrix priorvariance(d, d,1);
	//priorvariance.eye();
	priorvariance = Global::Psi*((Global::kappa + 1) / ((Global::kappa)*Global::eta));
	priormean = Global::mu0;
	

	IWishart w(priorvariance,d+2);
	w.rnd().print();
	w.rnd().print();
	w.rnd().print();
	vector<Normal> restaurant;
	ThreadPool tpool(thread::hardware_concurrency());
	Normal nr(priormean,priorvariance);
	Vector& r = nr.rnd();

	kmeans(ds);

	system("pause");
}
