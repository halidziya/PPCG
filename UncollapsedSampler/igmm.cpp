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
#include "Restaurant.h"
using namespace std;


int MAX_SWEEP = 500;
int BURNIN = 100;
int SAMPLE = (MAX_SWEEP - BURNIN) / 10; // Default value is 10 sample + 1 post burnin
char* result_dir = "";
int MAX_TABLE = 100;



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
		MAX_TABLE = atoi(argv[5]);
	if (argc > 6)
		result_dir = argv[6];
	//SAMPLE = (MAX_SWEEP - BURNIN) / 10; // Default value
	//if (argc>7)
	//	SAMPLE = atoi(argv[7]);

	step();
	
	string ss(result_dir);
	printf("size of int %d", sizeof(int));
	printf("Reading...\n");
	DataSet ds(datafile, priorfile, configfile);
	n = ds.data.r;
	d = ds.data.m;
	Vector priormean(d);
	Matrix priorvariance(d, d,1);
	//priorvariance.eye();
	priorvariance = Psi*((kappa + 1) / ((kappa)*T_eta));
	priormean = mu0;
	
	

	Vector labels = kmeans(ds);
	
	precomputeGammaLn(2 * n + 100*d);
	Stut stt(priormean, priorvariance, T_eta);
	loglik0 = stt.likelihood(ds.data);
	Restaurant r(ds, MAX_TABLE);
	r.createTables(labels);
	ThreadPool tpool(thread::hardware_concurrency());



	r.getInfo();
	for (auto i = 0; i < MAX_SWEEP; i++)
	{
		r.resetStats();
		for (auto i = 0; i < tpool.numthreads; i++) {
			tpool.submit(r);
		}
		tpool.waitAll();
		r.samplePosteriors();
		if (i % 100 == 0)
		{
			printf("\n\nITER  : %d\n\n", i);
			// r.getInfo();
			flush(cout);
		}
	}
	r.getInfo();

	string s(result_dir);
	ofstream restfile(s.append("_igmm.rest"), ios::out | ios::binary);

	restfile << r;
	restfile.close();
}
