#pragma once
#include "DataSet.h"
#include "Table.h"
#include "ThreadPool.h"
#include "Dirichlet.h"

extern Vector loglik0;


class Restaurant : public Task{
public:
	DataSet& ds;
	vector<list<Table>::iterator> labels;
	list<Table> tables;
	vector<list<int>> collapsedList;
	void run(int id);
	void resetStats();
	void samplePosteriors();
	void createTables(Vector& labels);
	Restaurant(DataSet& ds);
	Dirichlet tabledist;
	Vector sticks;
	int maxtable;
	void getInfo();
	friend ostream& operator<<(ostream& os,Restaurant& r);
	list<Table> sampleCollapsed(list<int> dataids);
	bool  collapsed = 1; //Last iteration is not collapsed
	atomic<int> taskid;
};
