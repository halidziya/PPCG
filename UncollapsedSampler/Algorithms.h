#include "Matrix.h"
#include "ThreadPool.h"
#include "DataSet.h"
class KmeansTask : public Task
{
public:
	Matrix data;
	Matrix means;
	Vector labels;
	Vector counts;
	Matrix sums;
	KmeansTask(Matrix& data, Matrix& means);
	void run(int id);
};

Vector kmeans(DataSet& ds,int k=3);