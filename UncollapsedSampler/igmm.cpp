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
int BURNIN = 300;
int STEP = (MAX_SWEEP - BURNIN) / 10; // Default value is 10 sample + 1 post burnin
char* result_dir = "";


double gammalnd(int x) // Actually works on x/2 
{
	double res = 0;
	for (auto i = 0; i < d; i++)  
	{
		res += gl_pc[x - i]; 
	}
	return (log(M_PI)*d*(d - 1) / 4) + res;
}



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
	if (argc > 5)
		result_dir = argv[5];
	//SAMPLE = (MAX_SWEEP - BURNIN) / 10; // Default value
	//if (argc>7)
	//	SAMPLE = atoi(argv[7]);

	step();
	
	string ss(result_dir);
	printf("Reading...\n");
	DataSet ds(datafile, priorfile, configfile);
	n = ds.data.r;
	d = ds.data.m;
	Vector priormean(d);
	Matrix priorvariance(d, d,1);
	//priorvariance.eye();
	priorvariance = Psi*((kappa + 1) / ((kappa)*T_eta));
	priormean = mu0;
	
	

	generator.seed(time(NULL));
	srand(time(NULL));


	Vector labels(n);// = kmeans(ds);
	labels.zero();
	double logpi = log(M_PI);
	precomputeGammaLn(2 * n + 100*d);


	Stut stt(priormean, priorvariance, T_eta);
	loglik0 = stt.likelihood(ds.data);
	Restaurant r(ds);
	r.createTables(labels);
	printf("\nInitial tables : %d\n", r.tables.size());
	ThreadPool tpool(thread::hardware_concurrency());

	int nlabelsample = ((MAX_SWEEP - BURNIN) / STEP);
	Matrix sampledLabels(nlabelsample,n);
	Vector likelihoods(MAX_SWEEP);
	for (auto iter = 0; iter < MAX_SWEEP; iter++)
	{

			r.resetStats();
			for (auto i = 0; i < tpool.numthreads; i++) {
				tpool.submit(r);
			}
			tpool.waitAll();
			r.samplePosteriors();
		//system("pause");
		
		if (iter % 20 == 0)
		{
			printf("\n\nITER  : %d\n\n", iter);
			printf("\nTables : %d\n", r.tables.size());
			// r.getInfo();
			flush(cout);
		}


		if (iter >= BURNIN && (iter - BURNIN)%STEP==0) {
			int li = (((iter - BURNIN) / STEP));
			for (auto i = 0; i < n;i++)
				sampledLabels(li)[i] = r.labels[i]->id;
		}

	double totallikelihood = 0;
		for (auto& table : r.tables) // Jchang's formula for joint marginal distribution
		{
			double s1 = kappa + table.totalpoints;
			Vector sum = table.sum[0];
			for (auto i = 1; i < tpool.numthreads; i++)
				sum = sum + table.sum[i];
			Vector& diff = (mu0 - (sum / table.totalpoints));
			Matrix& outer = (diff >> diff)*((table.totalpoints *kappa) / (s1));
			Matrix ss = Psi + outer;
			for (auto i = 0; i < tpool.numthreads; i++)
				ss = ss + table.scatter[i];
			ss = ss / (table.totalpoints + m);
			totallikelihood +=  -0.5*table.totalpoints*d*logpi - gammalnd(m) + gammalnd(table.totalpoints + m) - (0.5*(table.totalpoints + m))*(d*log(table.totalpoints + m)
				+ 2*ss.chol().sumlogdiag()) + (0.5*m)*(d*log(m) + 2*(Psi/m).chol().sumlogdiag()) - 0.5*d*log((table.totalpoints + kappa) / kappa) + log(gamma) + gl_pc[table.totalpoints*2];
		}
		likelihoods[iter] = totallikelihood + gl_pc[gamma*2] - gl_pc[(gamma+n)*2];
	}


	



	string s(result_dir);
	ofstream restfile(s.append("_igmm.rest"), ios::out | ios::binary);

	string s2(result_dir);
	ofstream labelfile(s2.append("_igmm.labels"), ios::out | ios::binary);

	string s3(result_dir);
	ofstream likefile(s3.append("_igmm.likelihood"), ios::out | ios::binary);

	restfile << r;
	restfile.close();

	labelfile << sampledLabels;
	labelfile.close();

	likefile << likelihoods;
	likefile.close();

}
