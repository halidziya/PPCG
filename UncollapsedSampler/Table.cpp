#include "Table.h"

double Table::loglikelihood(Vector & v)
{
		return datadist.likelihood(v);
}

void Table::calculatePosteriors()
{
	Vector ssum(d);
	Matrix sscatter(d,d);
	
	ssum.zero();
	sscatter.zero();

	calculateTotalPoints();
	
	if (totalpoints != 0) {
	for (auto i = 0; i < sum.size(); i++) {
		ssum = ssum + sum[i];
		sscatter = sscatter + scatter[i];
	}
		double divider = 1.0/(kappa+ totalpoints);
		Vector& mean = ssum / totalpoints;
		//sscatter = sscatter - (mean.outer(mean))*totalpoints;
		Vector& diff = (mean-mu0);
		posteriorcov = IWishart(Psi + sscatter + diff.outer(diff)*divider*kappa*totalpoints , m + totalpoints);
		Matrix& sigma = posteriorcov.rnd();
		
		posteriormean = Normal((ssum + mu0*kappa)*divider, sigma*divider);
		datadist = Normal(posteriormean.rnd(), sigma);

		/*
		Vector& mean = ssum / totalpoints;
		//sscatter = sscatter - (mean.outer(mean))*totalpoints;
		posteriorcov = IWishart(Psi + sscatter, m + totalpoints);
		Matrix& sigma = posteriorcov.rnd();
		Matrix& sigmainv = sigma.inverse();
		Matrix& newcov = (Psi.inverse() + sigmainv*totalpoints).inverse();
		posteriormean = Normal(newcov*(sigmainv*ssum + Psi.inverse()*mu0), newcov);
		datadist = Normal(posteriormean.rnd(), sigma);
		*/
	}
	//printf("--%d--\n", predictive);
}

void Table::resetStats()
{
	
	int nthd = thread::hardware_concurrency();
	sum.resize(nthd);
	scatter.resize(nthd);
	likelihood.resize(nthd);
	npoints.resize(nthd);
	for (auto i = 0; i < sum.size(); i++) {
		sum[i] = zeros(d);
		scatter[i] = zeros(d, d);
		likelihood[i] = 0;
		npoints[i] = 0;
	}
	totalpoints = 0;
}

void Table::addPoint(Vector & v,int id)
{
	this->sum[id] = this->sum[id] + v;
	this->npoints[id] += 1;
	Vector& diff = v - datadist.mu;
	this->scatter[id] = this->scatter[id] + diff.outer(diff);
}

//Used in collapsed version
void Table::cremovePoint(Vector & v)
{
	int threadid = 0;
	 // Get abstract of the output buffer
	//sampleScatter = (sampleScatter - (diff >> diff)*(npoints / (npoints - 1.0)));
	//sampleMean = ((sampleMean * (npoints / (npoints - 1.0))) - v *(1.0 / (npoints - 1.0)));
	this->sum[threadid] = this->sum[threadid] - v;
	Vector& sampleMean = sum[threadid] / npoints[threadid];
	Vector& diff = (v - sampleMean);
	this->scatter[threadid] = this->scatter[threadid] + (diff >> diff)*(npoints[threadid] / (npoints[threadid] - 1.0));
	this->npoints[threadid] -= 1;
}

//Used in collapsed version
void Table::caddPoint(Vector & v)
{
	int threadid = 0;
	this->npoints[threadid] += 1;
	Vector& sampleMean = sum[threadid] / (npoints[threadid]-1);
	if (npoints[threadid] == 1)
		sampleMean <= mu0;
	Vector& diff = (v - sampleMean);
	this->scatter[threadid] = this->scatter[threadid] + (diff >> diff)*((npoints[threadid] -1.0)/ npoints[threadid]);
	this->sum[threadid] = this->sum[threadid] + v;
}


void Table::calculateTotalPoints()
{
	totalpoints = 0;
	for (auto i = 0; i < sum.size(); i++)
		totalpoints = totalpoints + npoints[i];
}

Table::Table()
{
	posteriormean = Normal(d);
	posteriorcov = IWishart(Psi,m);
	datadist     = Normal(d); 
	resetStats();
}


ostream& operator<<(ostream& os, const Table& t)
{

	//printf("%d\n", t.totalpoints);
	os.write((char*) &t.totalpoints,sizeof(int));
	os << t.datadist.mu;
	os << t.datadist.cholsigma;
	return os;
}