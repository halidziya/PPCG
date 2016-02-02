#pragma once
#include "Matrix.h"
#include "Global.h"
#include "Stut.h"

class DataSet 
{
public:
	DataSet(char* datafile,char* priorfile,char* configfile);
	Matrix data;
	vector<Matrix> chunks;
	Matrix prior;
	~DataSet(void);
};

