/*
 *  Copyright 2012 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <math.h> 
#include <string.h>
//#include <openacc.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXPOINT 1000
#define MAXTABLE 100


int main(int argc, char** argv)
{
	int npoints = 1;
	int ntables = 1;
	int maxiter = 1;
	int d = 1;
	double* likemat= (double*) malloc(sizeof(double)*ntables*npoints);

	for (auto iter = 0; iter < maxiter; iter++)
	{
	#pragma acc kernels
	for (auto i = 0; i < npoints; i++)
	{
		for (auto j = 0; j < ntables; j++)
		{
			double likelihood = 0;
			for (auto di = 0; di < d; di++)
				likelihood += 0; // data - mu
			likemat[i+j*npoints] = likelihood;
		}
	}
	}
	system("pause");
}

