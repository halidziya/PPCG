#include "FastMat.h"

#include "util.h"
#include <iostream>
#include <fstream>
#include <list>
#include "Algorithms.h"
#include "Table.h"
#include "Restaurant.h"
using namespace std;

int MAX_SWEEP = 1500;
int BURNIN = 1000;
int SAMPLE = 20;
int STEP = (MAX_SWEEP - BURNIN) / SAMPLE;
char* result_dir = "./";

double gamlnd(int x) // Actually works on x/2 
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




	if (argc < 2) {
		printf("Usage: igmm.exe datafilename [priorfilename] [paramfilename] [MAX_SWEEP] [BURNIN] [result_dir] [SAMPLE]");
		exit(-1);
	}
	
	if (argc>4)
		MAX_SWEEP = atoi(argv[4]);
	if (argc>5)
		BURNIN = atoi(argv[5]);
	if (argc > 6)
		result_dir = argv[6];
	else
	{
		string str(argv[1]);
		result_dir = (char*)str.substr(0, str.find_last_of("/\\")).c_str(); // Datafile folder
	}
	if (argc > 7)
	{
		SAMPLE = atoi(argv[7]);
	}
	STEP = (MAX_SWEEP - BURNIN) / SAMPLE;
	if (BURNIN >= MAX_SWEEP | STEP == 0) // Housekeeping
	{
		BURNIN = MAX_SWEEP - 2;
		SAMPLE = 1; STEP = 1;
	}

	
	
	string ss(result_dir);
	printf("Reading...\n");
	nthd = thread::hardware_concurrency();
	DataSet ds(argc,argv);
	n = ds.data.r;
	d = ds.data.m;
	Vector priormean(d);
	Matrix priorvariance(d, d,1);
	T_eta = m - d + 1;
	//priorvariance.eye();
	priorvariance = Psi*((kappa + 1) / ((kappa)*T_eta));
	priormean = mu0;
	
	

	generator.seed(time(NULL));
	srand(time(NULL));


	Vector labels(n);// = kmeans(ds);
	labels.zero();
	double logpi = log(M_PI);
	precomputegamLn(2 * n + 100*d);


	Stut stt(priormean, priorvariance, T_eta);
	loglik0 = stt.likelihood(ds.data);
	Restaurant r(ds);
	r.createTables(labels);
	printf("\nInitial tables : %d\n", r.tables.size());
	ThreadPool tpool(thread::hardware_concurrency());
	Matrix     sampledLabels(SAMPLE, n);
	Vector likelihoods(MAX_SWEEP);
	PILL_DEBUG;
	step();
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


		//reid(r.tables);
		if (((MAX_SWEEP - iter - 1) % STEP) == 0 && iter >= BURNIN)
		{
			int sampleno = (MAX_SWEEP - iter - 1) / STEP;
			if (sampleno < SAMPLE)
			{
				for (auto i = 0; i < n; i++)
					sampledLabels(sampleno)[i] = r.labels[i]->id;
			}
		}

		//if (iter >= BURNIN && (iter - BURNIN)%STEP==0) {
		//	int li = (((iter - BURNIN) / STEP));
		//	for (auto i = 0; i < n;i++)
		//		sampledLabels(li)[i] = r.labels[i]->id;
		//}

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
			totallikelihood +=  -0.5*table.totalpoints*d*logpi - gamlnd(m) + gamlnd(table.totalpoints + m) - (0.5*(table.totalpoints + m))*(d*log(table.totalpoints + m)
				+ 2*ss.chol().sumlogdiag()) + (0.5*m)*(d*log(m) + 2*(Psi/m).chol().sumlogdiag()) - 0.5*d*log((table.totalpoints + kappa) / kappa) + log(gam) + gl_pc[table.totalpoints*2];
		}
		likelihoods[iter] = totallikelihood + gl_pc[gam*2] - gl_pc[(gam+n)*2];
	}


	


	step();
	
	string s(result_dir);
	ofstream restfile(s.append("_igmm.rest"), ios::out | ios::binary);

	//string s2(result_dir);
	//ofstream labelfile(s2.append("_igmm.labels"), ios::out | ios::binary);

	string s2(result_dir);
	ofstream labelfile(s2.append("Labels.matrix"), ios::out | ios::binary);

	string s3(result_dir);
	ofstream likefile(s3.append("_igmm.likelihood"), ios::out | ios::binary);

	restfile << r;
	restfile.close();

	labelfile << sampledLabels;
	labelfile.close();

	likefile << likelihoods;
	likefile.close();

}
