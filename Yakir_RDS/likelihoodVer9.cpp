#include <vector>
#include <limits>
#include <iostream>
#include <math.h>
#include <R.h>

using namespace std;

extern "C"{
	void likelihood(  int *sample,
					  int *Sij,
					  int *S,
					  double *c,
					  double *theta,
					  double *Nj,
					  double *constant,
					  int *observed_degrees,
					  int *n,
					  int *N,
					  int *N_observed,
					  bool *arc,
					  double *result ) {				
		
		bool			first = true;
		
		*result=0.0;
		int temp_degree=0;

		 //actual code
		for(int i=0; i < *n; ++i)  {
			for(int j=0; j < *N_observed; ++j) {
				temp_degree = observed_degrees[j];
				*result += (sample[i] == temp_degree ) ?
					(log(*c) + log((double)S[i]) + log(Nj[temp_degree-1] - Sij[j + i*(*N_observed) ]) + *theta * log((double)temp_degree)) :
					(log(1.0 - *c * S[i] * (Nj[temp_degree-1]-Sij[j + i*( *N_observed)] ) * pow(double(temp_degree), *theta)));
			}
			if(sample[i]>0 && first) 
				first = false;
			}
		
		if(isnan(*result))
			*result = -numeric_limits<double>::infinity();			
	}
}
