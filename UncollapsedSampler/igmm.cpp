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
#include "Table.h"
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
	
	

	Vector labels = kmeans(ds);
	Vector loglik0;
	precomputeGammaLn(2 * n + 10);
	Stut stt(priormean, priorvariance, ds.m);
	loglik0 = stt.likelihood(ds.data);

	//UncollapsedSampler(ds, labels);
	Normal nr(v({ 1.0,1.0,1.0 }), eye(d));
	cout << nr.likelihood(zeros(d));
	auto x = { 1,1,1 };

	


	system("pause");
}
