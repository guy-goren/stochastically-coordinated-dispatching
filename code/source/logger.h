#ifndef LOGGER_H_
#define LOGGER_H_

#include "defs.h"

class Logger
{

public:

	Logger(int simulation_time, int sample_resolution);
	~Logger();

	void updateMQL(int* queueLengths, int n, int time);
	void updateJCT(int job_creation_time, int time);

	void updateGDT(long long getDestinations_microsec_exectime);

	void printStatus(int time, int ps_resolution, int current_sum);

	void writeResultsToFiles(string fn, int time);

private:

	int st_;
	int sr_;

	int nsu_;
	
	// store historgam of job conpletion times (jct's - exact)
	unordered_map<int, int> jct_map_;

	// store the evolution of the time averaged queue lengths (sampled)
	vector<long double> ta_queue_lengths_;

	// cumulative sum of jct's (exact)
	long double job_jcts_cumsum_, completed_jobs_;

	// cumulative sum of queue lengths (exact)
	long double queue_lengths_cumsum_;

	// store historgam of getDestinations execution time
	unordered_map<long long, int> gdt_map_;

};

#endif

