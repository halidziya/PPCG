#pragma once
#include "Normal.h"
#include "IWishart.h"




class Table {
public:
	Normal posteriormean;
	IWishart posteriorcov;
	Normal datadist;
	double stick;

	vector<Vector> sum; // Sum for each thread
	vector<Matrix> scatter; // Scatter for each thread
	vector<double> likelihood; // Likelihood buffer for each thread, It avoids likelihood vector generation for each point
	vector<int>    npoints; // Number of data points for each thread
	int totalpoints;
	static double eta;
	static double kappa;
	double loglikelihood(Vector& v, double predictivelikelihood=-INFINITY);
	void calculatePosteriors();
	void resetStats();
	void addPoint(Vector& v, int threadid);
	bool predictive; // Use predictive distribtion
	Table();
	friend ostream& operator<<(ostream& os, const Table& t);
};