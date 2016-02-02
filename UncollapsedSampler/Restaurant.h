#pragma once
#include "DataSet.h"
#include "Table.h"
#include "ThreadPool.h"
#include "Dirichlet.h"
extern Vector loglik0;


class Restaurant : public Task{
public:
	DataSet& ds;
	Vector labels;
	vector<Table> tables;
	void run(int id);
	void resetStats();
	void samplePosteriors();
	void createTables(Vector& labels);
	Restaurant(DataSet& ds,int ntables);
	Dirichlet tabledist;
	Vector sticks;
	int maxtable;
	void getInfo();
	friend ostream& operator<<(ostream& os, const Restaurant& r);
};
