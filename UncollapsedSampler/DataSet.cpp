#include "DataSet.h"
#include <iostream>
#include <fstream>

using namespace std;

DataSet::DataSet(char* datafile,char* priorfile,char* configfile)
{
	Matrix conf;
	Matrix labels;
	data.readBin(datafile);
	prior.readBin(priorfile);
	conf.readBin(configfile);
	n = data.r;
	d =  conf(0)(0);
	m =  conf(0)[1];
	kappa  = conf(0)[2];
	gamma  = conf(0)[3];

	int nthd = thread::hardware_concurrency();
	init_buffer(nthd, d);
	printf("Threads %d\n", nthd);

	Psi = prior;
	Psi.r = d; // Last row is shadowed
	Psi.n = d*d;
	mu0 = prior(d).copy();
	T_eta = m - d + 1;

	int chunksize = ceil((double)data.r / nthd);
	chunks.resize(nthd);
	for (auto i = 0; i < nthd; i++)
	{
		chunks[i] <= data.submat(i*chunksize,(i+1)*chunksize,0,d);
	}

}


DataSet::~DataSet(void)
{
	Psi.r = d+1; // Last row is shadowed
	Psi.n = d*(d+1);
}


