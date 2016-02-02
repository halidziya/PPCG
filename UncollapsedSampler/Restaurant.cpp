#include "Restaurant.h"


Vector loglik0;
void Restaurant::run(int id)
{
	
	// Use thread specific buffer
	matbuffer.threadid = id;
	buffer.threadid = id;
	absbuffer.threadid = id;

	id = id - 1; // Oth thread is the main process
	
	Matrix& chunk = ds.chunks[id];
	double maxlikelihood = 0, sum=0,val = 0,i,j;
	int chunksize = ds.chunks[0].r;
	//printf("id : %d  Last one : %d , chunk size : %d datasize : %d \n",id, id*chunksize + data.r - 1, chunksize, data.r);
	for ( i = 0; i < chunk.r; i++)
	{
		Vector& data = chunk(i);
		for ( j = 0; j < tables.size(); j++)
		{
			//printf("%.4f\n", tables[j].loglikelihood(data(i),-5));
			tables[j].likelihood[id] = tables[j].loglikelihood(data,loglik0[id*chunksize + i])+log(sticks[j]); //Log likelihoods sum
		}
		//printf("\n==\n");

		maxlikelihood = tables[0].likelihood[id];
		for (auto& table : tables)
		{
			if (table.likelihood[id] > maxlikelihood)
				maxlikelihood = table.likelihood[id];
		}

		sum = 0;
		for (auto& table : tables)
		{
				table.likelihood[id] = exp(table.likelihood[id] - maxlikelihood);
				sum += table.likelihood[id];
		}

		val = urand()*sum;
		for (j = 0; j < tables.size(); j++)
		{
			if (tables[j].likelihood[id] >= val)
			{
				break; // Find it
			}
			val -= (tables[j].likelihood[id]);
		}
		
		tables[j].addPoint(data, id);

		labels[id*chunksize +  i] = j;
	}

}

void Restaurant::resetStats()
{
	for (auto& table : tables)
		table.resetStats();

}

void Restaurant::samplePosteriors()
{
	Vector alpha(maxtable, 1, gamma);
	for (auto i = 0; i < maxtable; i++)
	{
		tables[i].calculatePosteriors();
		alpha[i] += tables[i].totalpoints;
	}
	tabledist = Dirichlet(alpha);
	sticks = tabledist.rnd();
}


void Restaurant::createTables(Vector& labels)
{
	srand(time(NULL));
	Vector& tableids = labels.unique();
	Vector alpha(maxtable, 1, gamma);
	tables.resize(maxtable);
	for (auto& table : tables)
		table.resetStats();
	for (auto i = 0; i < labels.n; i++)
		tables[labels[i]].addPoint(ds.data(i),0);
	for (auto i = 0; i < maxtable; i++)
	{
		tables[i].calculatePosteriors();
		alpha[i] += tables[i].totalpoints;
	}
	tabledist = Dirichlet(alpha);
	printf("Alpha sum %d", alpha.sum());
	sticks = tabledist.rnd();
}


Restaurant::Restaurant(DataSet& ds,int ntables) : ds(ds) {
	maxtable = ntables;
	labels.resize(ds.data.r);
}

void Restaurant::getInfo()
{
	printf("+++++++++++RESTAURANT++++++++++\n\n");
	int tableid = 0;
	for (auto& table : tables)
	{
		table.datadist.mu.print();
		for (auto& i : table.npoints)
			printf("%d ", i);
		printf("==> %d %.3f\n", table.totalpoints, sticks[tableid]);
		tableid++;
	}
}

ostream & operator<<(ostream & os, const Restaurant & r)
{
	int a = r.tables.size();
	os.write((char*)&a, sizeof(int));
	for (const Table& table : r.tables)
		os << table;
	os << r.labels;
	return os;
}
