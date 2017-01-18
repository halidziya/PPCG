#pragma once
#include "FastMat.h"
#include "GMMBase.h"



class Table {
public:
	Normal posteriormean;
	IWishart posteriorcov;
	Normal datadist;

	vector<Vector> sum; // Sum for each thread
	vector<Matrix> scatter; // Scatter for each thread
	vector<double> likelihood; // Likelihood buffer for each thread, It avoids likelihood vector generation for each point
	vector<int>    npoints; // Number of data points for each thread
	int totalpoints;
	double loglikelihood(Vector& v);
	void calculatePosteriors();
	void resetStats();
	void addPoint(Vector& v, int threadid);
	void cremovePoint(Vector& v); // Collapsed Table
	void caddPoint(Vector& v); // Collapsed Table
	void calculateTotalPoints();
	Table();
	friend ostream& operator<<(ostream& os, Table& t);
	int id; // Used in collapsed part, uncollapsed one does not need for now
};