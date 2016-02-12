#pragma once
#include "DataSet.h"
#include "Table.h"
#include "ThreadPool.h"
#include "Dirichlet.h"
#include <atomic>
extern Vector loglik0;


class Restaurant : public Task{
public:
	DataSet& ds;
	Vector labels;
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
	boolean  collapsed = 1; //Last iteration is not collapsed
	atomic<int> taskid;
};
