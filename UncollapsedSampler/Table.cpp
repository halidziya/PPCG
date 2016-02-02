#include "Table.h"

double Table::eta;
double Table::kappa;



double Table::loglikelihood(Vector & v,double predictivelikelihood)
{
	if (!predictive) // Occupied Table
	{
		return datadist.likelihood(v);
	}
	else
		return predictivelikelihood;
}

void Table::calculatePosteriors()
{
	Vector ssum(d);
	Matrix sscatter(d,d);
	totalpoints=0;
	ssum.zero();
	sscatter.zero();
	for (auto i = 0; i < sum.size(); i++)
		totalpoints = totalpoints + npoints[i];
	
	if (totalpoints != 0) {
	for (auto i = 0; i < sum.size(); i++) {
		ssum = ssum + sum[i];
		sscatter = sscatter + scatter[i];
	}

		predictive = false;
		Vector& mean = ssum / totalpoints;
		//sscatter = sscatter - (mean.outer(mean))*totalpoints;
		posteriorcov = IWishart(Psi + sscatter, m + totalpoints);
		Matrix& sigma = posteriorcov.rnd();
		Matrix& sigmainv = sigma.inverse();
		Matrix& newcov = (Psi.inverse() + sigmainv*totalpoints).inverse();
		posteriormean = Normal(newcov*(sigmainv*ssum + Psi.inverse()*mu0), newcov);
		datadist = Normal(posteriormean.rnd(), sigma);
	}
	else
		predictive = true; // access to student T likelihood
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

}

void Table::addPoint(Vector & v,int id)
{
	this->sum[id] = this->sum[id] + v;
	this->npoints[id] += 1;
	Vector& diff = v - datadist.mu;
	this->scatter[id] = this->scatter[id] + diff.outer(diff);
}

Table::Table()
{
	posteriormean = Normal(d);
	posteriorcov = IWishart(Psi,eta);
	datadist     = Normal(d); 
	resetStats();
}


ostream& operator<<(ostream& os, const Table& t)
{

	printf("%d\n", t.totalpoints);
	os.write((char*) &t.totalpoints,sizeof(int));
	os << t.datadist.mu;
	os << t.datadist.cholsigma;
	return os;
}