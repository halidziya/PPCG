#include "Restaurant.h"


Vector loglik0;
void Restaurant::run(int id)
{
	
	// Use thread specific buffer
	matbuffer.threadid = id;
	buffer.threadid = id;
	absbuffer.threadid = id;

	int taskid = this->taskid++; // Oth thread is the main process
	
	Matrix& chunk = ds.chunks[taskid];
	double maxlikelihood = 0, sum = 0, val = 0, newtablelikelihood = 0;
	int i, j,dataid;
	int chunksize = ds.chunks[0].r;
	
	//printf("id : %d  Last one : %d , chunk size : %d datasize : %d \n",id, id*chunksize + data.r - 1, chunksize, data.r);
	for ( i = 0; i < chunk.r; i++)
	{
		dataid = taskid*chunksize + i;
		Vector& data = chunk(i);
		j = 0;
		for (auto& table : tables)
		{
			//printf("%.4f\n", tables[j].loglikelihood(data(i),-5));
			table.likelihood[taskid] = table.loglikelihood(data)+log(sticks[j]); //Log likelihoods sum
			j++;
		}
		newtablelikelihood = loglik0[dataid] + log(sticks[j]);
		//printf("\n==\n");

		maxlikelihood = tables.begin()->likelihood[taskid];
		for (auto& table : tables)
		{
			if (table.likelihood[taskid] > maxlikelihood)
				maxlikelihood = table.likelihood[taskid];
		}
		


		sum = 0;
		//printf("Point %d Label : %d ",dataid,labels[dataid]);
		//data.print();
		for (auto& table : tables)
		{
				//printf("%.3f\t", table.likelihood[id]);
				table.likelihood[taskid] = exp(table.likelihood[taskid] - maxlikelihood);
				sum += table.likelihood[taskid];
		}
		//printf("%.3f\t", newtablelikelihood);
		//printf("\n\n");

		if (collapsed) {
			// New table
			newtablelikelihood = newtablelikelihood - maxlikelihood;
			sum += exp(newtablelikelihood);
		}

		val = urand()*sum;
		j = 0;
		list<Table>::iterator lasttable=tables.end();
		for (auto table = tables.begin(); table != tables.end();table++)
		{
			if (table->likelihood[taskid] >= val)
			{
				lasttable = table;
				break; // Find it
			}
			val -= (table->likelihood[taskid]);
			j++;
		}

		if (collapsed && (lasttable == tables.end()))
		{
			labels[dataid] = -1; // Not yet assigned
			collapsedList[taskid].push_back(dataid);
		}
		else
		{
			lasttable->addPoint(data, taskid);
			labels[dataid] = j;
			//printf("Point : %d in %d with %d -> %d\n",dataid, taskid,i, j);
		}
	}

}

void Restaurant::resetStats()
{
	taskid = 0;
	for (auto& table : tables)
		table.resetStats();
	for (auto& alist : collapsedList)
		alist.clear();
}

void Restaurant::samplePosteriors()
{
	for (auto& t : tables)
	{
		t.calculateTotalPoints();
		//printf("\n|%d| ", t.totalpoints);
	}




	if (collapsed) {
		for (auto iter = tables.begin(); iter != tables.end();)
			if (iter->totalpoints == 0)
			{
				iter = tables.erase(iter);
			}
			else
				iter++;

		list<int> collection;
		for (auto& l : collapsedList)
			collection.insert(collection.end(), l.begin(), l.end());

		if (collection.size() > 0) { // New tables created
			list<Table> newtables = sampleCollapsed(collection);
			tables.insert(tables.end(), newtables.begin(), newtables.end());
		}
		
	}


	Vector alpha(tables.size() + 1, 1, 0);
	int  i = 0;
	for (auto& t : tables)
	{
		t.calculatePosteriors();
		alpha[i] += t.totalpoints;
		i++;
	}
	alpha[i] = gamma; // Non parametric part

	tabledist = Dirichlet(alpha);
	sticks = tabledist.rnd();
}


void Restaurant::createTables(Vector& labels)
{
	srand(time(NULL));
	Vector& tableids = labels.unique();
	Vector alpha(tableids.n+1, 1, 0);
	tables.resize(tableids.n);
	for (auto& table : tables)
		table.resetStats();

	vector<list<Table>::iterator> tablelist;
	for (auto it = tables.begin(); it != tables.end(); it++)
		tablelist.push_back(it);

	for (auto i = 0; i < labels.n; i++) {
		tablelist[labels[i]]->caddPoint(ds.data(i));
		this->labels[i] = labels[i];
	}
	
	int  i = 0;
	for (auto& t : tables)
	{
		t.calculatePosteriors();
		alpha[i] += t.totalpoints;
		i++;
	}
	alpha[i] = gamma;
	collapsedList.resize(thread::hardware_concurrency());
	tabledist = Dirichlet(alpha);
	printf("Alpha sum %d", alpha.sum());
	sticks = tabledist.rnd();
}


Restaurant::Restaurant(DataSet& ds) : ds(ds) {
	labels.resize(ds.data.r);
	collapsed = true;
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

list<Table> Restaurant::sampleCollapsed(list<int> dataids)
{
	list<Table> ctables(1); //First table
	double sum, maxlikelihood,val;
	int i, j;
	//for (auto iter = 0; iter < 5; iter++) // 5 inner collapsed iter
	//{
	//
	//}
	int newtableid = tables.size();
	list<Table>::iterator it = ctables.begin();
	list<list<Table>::iterator> collapsedLabels;
	it->resetStats();
	// Initially all in same table
	for (auto dataid : dataids) {
		
		it->caddPoint(ds.data(dataid)); // Only 1 thread
		labels[dataid] = newtableid;
	}
	
	for (auto iter = 0; iter < 1;iter++){
	if (dataids.size()>1)
	for (auto dataid : dataids) {
		/*
		Vector& v = ds.data(dataid);
		it->cremovePoint(v);
		
		if (it->npoints[0] == 0)
			ctables.erase(it);
		/*	void calculateDist()
	{
		double s1 = Global::kappa + npoints;
		dist.eta = m + 1 + this->npoints - d;
		dist.mu = (sampleMean*(npoints / s1) + mu0 * (kappa / s1));
		Vector& diff = mu0 - sampleMean;
		dist.cholsigma =
			((Psi + sampleScatter +
				((diff) >> (diff))
				*(npoints*kappa / (s1)))
				*((s1 + 1) / (s1*dist.eta))).chol();
		dist.calculateNormalizer();
		}
	*/
	/*
		for (auto& table : ctables)
		{
			int npoints = table.npoints[0];
			double s1 = kappa + npoints;
			Vector& diff = (mu0 - table.sum[0] / npoints);
			Matrix& outer = (diff >> diff)*((npoints *kappa)/(s1));
			table.likelihood[0] = table.npoints[0] * Stut((table.sum[0] + mu0*kappa) / (s1), (Psi + table.scatter[0] + outer) * ((s1 + 1) / (s1*(m+npoints-d+1))), m + npoints - d + 1).likelihood(v);
		}
		double newtable = gamma*loglik0[dataid];
		

		maxlikelihood = newtable;
		for (auto& table : ctables)
		{
			if (table.likelihood[0] > maxlikelihood)
				maxlikelihood = table.likelihood[0];
		}

		sum = 0;
		for (auto& table : ctables)
		{
			table.likelihood[0] = exp(table.likelihood[0] - maxlikelihood);
			sum += table.likelihood[0];
		}

		newtable = newtable - maxlikelihood;
		sum += exp(newtable);
		val = urand()*sum;
		j = 0;
		list<Table>::iterator lasttable = ctables.end();
		for (auto table = ctables.begin(); table != ctables.end(); table++)
		{
			if (table->likelihood[0] >= val)
			{
				lasttable = table;
				break; // Find it
			}
			val -= (table->likelihood[0]);
			j++;
		}

		if (lasttable == ctables.end())
		{
			Table t;
			t.resetStats();
			t.caddPoint(v);
			ctables.push_back(t);
			labels[dataid] = -1;
		}
		else
		{
			lasttable->caddPoint(v);
			labels[dataid] = -1;
		}
		*/
	}

	/*for (auto iter = ctables.begin(); iter != ctables.end();)
		if (iter->npoints[0] == 0)
		{
			iter = ctables.erase(iter);
		}
		else
			iter++;*/
			
	}
	
	//printf("==> %d", newtableid);

	return ctables;
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
