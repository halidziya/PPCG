#include "Algorithms.h"


Vector kmeans(DataSet& ds,int k)
{
	Matrix means(k, ds.d);
	Vector counts(k);
	Vector results(ds.n);
	ThreadPool tpool(thread::hardware_concurrency());
	vector<KmeansTask> ktasks;
	for (auto iter = 0; iter < 100; iter++){
		
		
		for (auto i = 0; i < k; i++) {
			means(i) = ds.data(i*(ds.n/k));
		}

		for (auto i = 0; i < tpool.numthreads; i++) {
			ktasks.push_back(KmeansTask(ds.chunks[i], means));
		}
		for (auto i = 0; i < tpool.numthreads; i++) {
			tpool.submit(ktasks[i]);
		}
		tpool.waitAll();

		means.zero();
		counts.zero();
		
		for (auto i = 0; i < tpool.numthreads; i++) {
			means = means + ktasks[i].sums;
			counts = counts + ktasks[i].counts;
		}

		for (auto i = 0; i < k; i++) {
			means(i) = means(i) / counts(i);
		}
	}

	int position = 0;
	for (auto i = 0; i < tpool.numthreads; i++) {
		results.put(position, ktasks[i].labels);
		position = position + ktasks[i].labels.n;
	}
	return results;
}

KmeansTask::KmeansTask(Matrix& data, Matrix& centers) :  data(data), means(centers)
{

}

void KmeansTask::run(int id)
{
	matbuffer.threadid = id;
	buffer.threadid = id;
	absbuffer.threadid = id;
	int n = data.r;
	int l;
	double mindist, dist;
	int k = means.r;
	int d = data.m;
	labels = Vector(n);
	counts = Vector(k);
	sums = Matrix(k,d);
	sums.zero();
	counts.zero();
	for (auto i = 0; i < n; i++)
	{
		mindist = my_infinity();
		for (auto j = 0; j < k; j++)
		{
			dist = (data(i) - means(j)).norm();
			if (dist < mindist) {
				mindist = dist;
				l = j;
			}
		}
		labels[i] = l;
		counts[l]++;
		sums(l) =  sums(l) +  data(i);
	}
}
