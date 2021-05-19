#ifndef RP_H_
#define RP_H_

#include "defs.h"

class RandomProcesses
{

public:

	RandomProcesses(int seed, int m, double* arrival_parameters, int n, double* departure_parameters);
	~RandomProcesses();

	int arrivals(int i);
	int departures(int i);

private:

	// servers - length of distribution vectors
	int n_;

	// dispatchers - length of distribution vectors
	int m_;

	// seed
	default_random_engine generator_;

	// arrivals
	poisson_distribution<int>** p_distribution_;

	// departures
	geometric_distribution<int>** g_distribution_;

};



#endif
