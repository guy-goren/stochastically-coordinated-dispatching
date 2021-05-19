#include "defs.h"
#include "simulations.h"
#include "dispatcher.h"

int main(int argc, char* argv[])
{
	
	// number of servers
	int n = atoi(argv[1]);

	// number of dispatchers
	int m = atoi(argv[2]);
	
	// simulation time-slots
	int time = atoi(argv[3]);

	// arrival parameters
	double* arrival_parameters = new double[m]();

	for (int i = 0; i < m; ++i)
	{
		arrival_parameters[i] = atof(argv[4 + i]);
	}

	// departure parameters
	double* departure_parameters = new double[n]();

	for (int i = 0; i < n; ++i)
	{
		departure_parameters[i] = atof(argv[4 + m + i]);
	}

	// simulation type
	string simulationType = argv[4 + m + n];

	// default parameters
	int d = 2;

	// policies
	if ((simulationType.compare("complete_information") != 0) && (simulationType.compare("wr") != 0) && (simulationType.compare("jiq") != 0) && (simulationType.compare("lsq") != 0))
	{
		cout << "received unsupported algorithm " << simulationType << endl;
		exit(1);
	}
	else if (simulationType.compare("lsq") == 0)
	{
		d = atoi(argv[4 + m + n + 1]);

		// the number of parameters depends on the algorithm
		assert(argc == (5 + n + m) && "wrong number of input parameters");
	}
	else
	{
		// the number of parameters depends on the algorithm
		assert(argc == (4 + n + m) && "wrong number of input parameters");
	}
	
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////

	int seed = 42;
	int fn_resolution = 100;

	/////////////////////////////////////////
	/////////////////////////////////////////

	// complete information
	if (simulationType.compare("complete_information") == 0)
	{	
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_htwf_slow");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_htwf_fast");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_htwf");
		complete_information_TWF(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_twf");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_jsq");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_hjsq");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_jfsq");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_hjfsq");
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_hjfsq_fast");
	}

	/////////////////////////////////////////
	/////////////////////////////////////////

	// weigted_random
	else if (simulationType.compare("wr") == 0)
	{
		complete_information(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_weighted_random");
	}

	/////////////////////////////////////////
	/////////////////////////////////////////

	// classical jiq
	else if (simulationType.compare("jiq") == 0)
	{
		classical_jiq(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_jiq");
		classical_jiq(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, "splittable_hjfiq");
	}

	/////////////////////////////////////////
	/////////////////////////////////////////

	// lsq
	else if (simulationType.compare("lsq") == 0)
	{
		lsq_sample(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, d, "splittable_jsq");
		lsq_sample_het(n, m, time, arrival_parameters, departure_parameters, seed, fn_resolution, d, "splittable_hjfsq");
	}

	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
	/////////////////////////////////////////

	delete[] arrival_parameters;
	delete[] departure_parameters;

	return 0;
	
}