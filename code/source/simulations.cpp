#include <iomanip>
#include <sstream>

#include "defs.h"
#include "logger.h"
#include "RandomProcesses.h"
#include "dispatcher.h"
#include "server.h"

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

string returnLoadString(const double* arrival_parameters, int m, const double* departure_parameters, int n)
{
	double sumOfArrivalRates = 0;
	for (int i = 0; i < m; ++i)
	{
		sumOfArrivalRates += arrival_parameters[i];
	}

	double sumOfDepartureRates = 0;
	for (int i = 0; i < n; ++i)
	{
		sumOfDepartureRates += departure_parameters[i];
	}

	double load = sumOfArrivalRates / sumOfDepartureRates;

	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << load;
	string l = stream.str();

	return l;
}

string returnProbabilityString(double probability)
{
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << probability;
	string p = stream.str();

	return p;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void complete_information(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType)
{

	// result file names
	string fn;
	fn.append(dispatcherType);
	fn.append("_");
	fn.append(to_string(n));
	fn.append("_ser_");
	fn.append(to_string(m));
	fn.append("_dis_");
	fn.append(returnLoadString(arrival_parameters, m, departure_parameters, n));
	fn.append("_load_");
	fn.append(to_string(seed));
	fn.append("_seed_");

	// we use randomization here
	default_random_engine local_generator(seed);

	// collect results
	Logger logger(time, fn_resolution);

	// arrivals and departures
	RandomProcesses rs(seed, m, arrival_parameters, n, departure_parameters);

	// servers initialization
	Server** servers = new Server*[n];
	for (int i = 0; i < n; i++)
	{
		servers[i] = new Server(i);
	}

	int* jobsQueueSizes;
	jobsQueueSizes = new int[n]();

	// Dispatchers initialization
	Dispatcher** dispatchers = new Dispatcher*[m];

	// instantiate according to dispatcherType
	if (dispatcherType.compare("splittable_htwf") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerState(i, m, n, departure_parameters, seed);
		}
	} 
	else if (dispatcherType.compare("splittable_htwf_slow") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new SlowHeterogenousDispatcherSplittableFullServerState(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_htwf_fast") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new FastHeterogenousDispatcherSplittableFullServerState(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_jsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerStateJSQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_jfsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerStateJFSQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_hjsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerStateHJSQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_hjfsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerStateHJFSQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_hjfsq_fast") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableFullServerStateFastHJFSQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_weighted_random") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableWeightedRandom(i, m, n, departure_parameters, seed);
		}
	}
	else
	{
		cout << "Dispatcher type received: " << dispatcherType << endl;
		assert("Unknown Dispatcher Type");
	}

	// Main loop
	for (int t = 1; t < time + 1; ++t)
	{

		// collect queue sizes from servers
		for (int srv = 0; srv < n; ++srv)
		{
			jobsQueueSizes[srv] = servers[srv]->getQueueSize();
		}

		// logger: print status every 1% of time
		int curr_sum = 0; for (int i = 0; i < n; i++) curr_sum += jobsQueueSizes[i];
		logger.printStatus(t, time / 100, curr_sum);

		// logger: collect qls
		logger.updateMQL(jobsQueueSizes, n, t);

		// arrivals
		for (int disp = 0; disp < m; ++disp)
		{

			auto t1 = std::chrono::high_resolution_clock::now();
			vector<DestinationJobsPairs> djPairs = dispatchers[disp]->getDestinations(jobsQueueSizes, rs.arrivals(disp));
			auto t2 = std::chrono::high_resolution_clock::now();

			auto getDestinations_nanosec_exectime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
			logger.updateGDT((long long)getDestinations_nanosec_exectime);

			for (int j = 0; j < (int)djPairs.size(); ++j)
			{
				int destination = djPairs[j].first;
				int destinationJobs = djPairs[j].second;

				// send jobs to destination
				servers[destination]->insertJobs(destinationJobs, t);
			}
		}

		// departures
		for (int i = 0; i < n; ++i)
		{
			int departures = rs.departures(i);
			vector<int> jobCreationTimes = servers[i]->serveJobsAndReturnCreationTimes(departures);
			for (int j = 0; j < (int)jobCreationTimes.size(); ++j)
			{
				// logger: collect jct
				logger.updateJCT(jobCreationTimes[j], t);
			}
		}

	}

	// logger: write results to a file
	logger.writeResultsToFiles(fn, time);

	// release memory

	delete[] jobsQueueSizes;

	for (int i = 0; i < n; i++)
	{
		delete servers[i];
	}
	delete[] servers;

	for (int i = 0; i < m; i++)
	{
		delete dispatchers[i];
	}
	delete[] dispatchers;

}

void complete_information_TWF(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType)
{

	// result file names
	string fn;
	fn.append(dispatcherType);
	fn.append("_");
	fn.append(to_string(n));
	fn.append("_ser_");
	fn.append(to_string(m));
	fn.append("_dis_");
	fn.append(returnLoadString(arrival_parameters, m, departure_parameters, n));
	fn.append("_load_");
	fn.append(to_string(seed));
	fn.append("_seed_");

	// we use randomization here
	default_random_engine local_generator(seed);

	// collect results
	Logger logger(time, fn_resolution);

	// arrivals and departures
	RandomProcesses rs(seed, m, arrival_parameters, n, departure_parameters);

	// servers initialization
	Server** servers = new Server * [n];
	for (int i = 0; i < n; i++)
	{
		servers[i] = new Server(i);
	}

	int* jobsQueueSizes;
	jobsQueueSizes = new int[n]();

	// Dispatchers initialization
	DispatcherTWF** dispatchers = new DispatcherTWF * [m];

	// instantiate according to dispatcherType
	if (dispatcherType.compare("splittable_twf") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherTWFSplittableFullServerState(i, m, n, departure_parameters, seed);
		}
	}
	else
	{
		cout << "Dispatcher type received: " << dispatcherType << endl;
		assert("Unknown Dispatcher Type");
	}

	// Main loop
	for (int t = 1; t < time + 1; ++t)
	{

		// collect queue sizes from servers
		for (int srv = 0; srv < n; ++srv)
		{
			jobsQueueSizes[srv] = servers[srv]->getQueueSize();
		}

		// logger: print status every 1% of time
		int curr_sum = 0; for (int i = 0; i < n; i++) curr_sum += jobsQueueSizes[i];
		logger.printStatus(t, time / 100, curr_sum);

		// logger: collect qls
		logger.updateMQL(jobsQueueSizes, n, t);

		// arrivals
		for (int disp = 0; disp < m; ++disp)
		{

			auto t1 = std::chrono::high_resolution_clock::now();
			vector<DestinationJobsPairs> djPairs = dispatchers[disp]->getDestinations(jobsQueueSizes, rs.arrivals(disp));
			auto t2 = std::chrono::high_resolution_clock::now();

			auto getDestinations_nanosec_exectime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
			logger.updateGDT((long long)getDestinations_nanosec_exectime);
			
			for (int j = 0; j < (int)djPairs.size(); ++j)
			{
				int destination = djPairs[j].first;
				int destinationJobs = djPairs[j].second;

				// send jobs to destination
				servers[destination]->insertJobs(destinationJobs, t);
			}
		}

		// departures
		for (int i = 0; i < n; ++i)
		{
			int departures = rs.departures(i);
			vector<int> jobCreationTimes = servers[i]->serveJobsAndReturnCreationTimes(departures);
			for (int j = 0; j < (int)jobCreationTimes.size(); ++j)
			{
				// logger: collect jct
				logger.updateJCT(jobCreationTimes[j], t);
			}
		}

	}

	// logger: write results to a file
	logger.writeResultsToFiles(fn, time);

	// release memory

	delete[] jobsQueueSizes;

	for (int i = 0; i < n; i++)
	{
		delete servers[i];
	}
	delete[] servers;

	for (int i = 0; i < m; i++)
	{
		delete dispatchers[i];
	}
	delete[] dispatchers;

}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void classical_jiq(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType)
{

	// result file names
	string fn;
	fn.append(dispatcherType);
	fn.append("_");
	fn.append(to_string(n));
	fn.append("_ser_");
	fn.append(to_string(m));
	fn.append("_dis_");
	fn.append(returnLoadString(arrival_parameters, m, departure_parameters, n));
	fn.append("_load_");
	fn.append(to_string(seed));
	fn.append("_seed_");

	// we use randomization here
	default_random_engine local_generator(seed);

	// for power of d - server randomization
	uniform_int_distribution<int> unifDistribution(0, m - 1);

	// collect results
	Logger logger(time, fn_resolution);

	// arrivals and departures
	RandomProcesses rs(seed, m, arrival_parameters, n, departure_parameters);

	// servers initialization
	ServerJIQ** servers = new ServerJIQ*[n];
	for (int i = 0; i < n; i++)
	{
		servers[i] = new ServerJIQ(i);
	}

	int* jobsQueueSizes;
	jobsQueueSizes = new int[n]();

	bool* QueueIdleTokenFree;
	QueueIdleTokenFree = new bool[n];
	for (int i = 0; i < n; i++)
	{
		QueueIdleTokenFree[i] = true;
	}

	// Dispatchers initialization
	DispatcherJIQ** dispatchers = new DispatcherJIQ *[m];

	// instantiate according to dispatcherType
	if (dispatcherType.compare("splittable_jiq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableJIQ(i, m, n, departure_parameters, seed);
		}
	}
	else if (dispatcherType.compare("splittable_hjfiq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableHetJFIQ(i, m, n, departure_parameters, seed);
		}
	}
	else
	{
		cout << "Dispatcher type received: " << dispatcherType << endl;
		assert("Unknown Dispatcher Type");
	}

	// Main loop
	for (int t = 1; t < time + 1; ++t)
	{

		// collect queue sizes from servers
		for (int srv = 0; srv < n; ++srv)
		{
			jobsQueueSizes[srv] = servers[srv]->getQueueSize();
		}

		// logger: print status every 1% of time
		int curr_sum = 0; for (int i = 0; i < n; i++) curr_sum += jobsQueueSizes[i];
		logger.printStatus(t, time / 100, curr_sum);

		// logger: collect qls
		logger.updateMQL(jobsQueueSizes, n, t);

		// arrivals
		for (int disp = 0; disp < m; ++disp)
		{
			vector<DestinationJobsPairs> djPairs = dispatchers[disp]->getDestinations(rs.arrivals(disp));
			for (int j = 0; j < (int)djPairs.size(); ++j)
			{
				int destination = djPairs[j].first;
				int destinationJobs = djPairs[j].second;

				// send jobs to destination
				servers[destination]->insertJobs(destinationJobs, t);

				// make sure idle-tokens is released
				QueueIdleTokenFree[destination] = true;
			}
		}

		// departures
		for (int i = 0; i < n; ++i)
		{
			int departures = rs.departures(i);
			vector<int> jobCreationTimes = servers[i]->serveJobsAndReturnCreationTimesJIQ(departures);
			for (int j = 0; j < (int)jobCreationTimes.size(); ++j)
			{
				// logger: collect jct
				logger.updateJCT(jobCreationTimes[j], t);
			}

			// update a dispatcher?
			if ((servers[i]->gotIdle()) && (QueueIdleTokenFree[i]))
			{
				int random_dispatcher = unifDistribution(local_generator);
				dispatchers[random_dispatcher]->updateIdleQueue(i);
				QueueIdleTokenFree[i] = false;
			}
		}

	}

	// logger: write results to a file
	logger.writeResultsToFiles(fn, time);

	// release memory

	delete[] QueueIdleTokenFree;
	delete[] jobsQueueSizes;

	for (int i = 0; i < n; i++)
	{
		delete servers[i];
	}
	delete[] servers;

	for (int i = 0; i < m; i++)
	{
		delete dispatchers[i];
	}
	delete[] dispatchers;

}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void lsq_sample(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, int d, string dispatcherType)
{

	// result file names
	string fn;
	fn.append(dispatcherType);
	fn.append("_lsq_");
	fn.append(to_string(d));
	fn.append("_chs_");
	fn.append(to_string(n));
	fn.append("_ser_");
	fn.append(to_string(m));
	fn.append("_dis_");
	fn.append(returnLoadString(arrival_parameters, m, departure_parameters, n));
	fn.append("_load_");
	fn.append(to_string(seed));
	fn.append("_seed_");

	// we use randomization here
	default_random_engine local_generator(seed);

	// for power of d - server randomization
	uniform_int_distribution<int> unifDistribution(0, n - 1);

	// collect results
	Logger logger(time, fn_resolution);

	// arrivals and departures
	RandomProcesses rs(seed, m, arrival_parameters, n, departure_parameters);

	// servers initialization
	ServerBasic** servers = new ServerBasic * [n];
	for (int i = 0; i < n; i++)
	{
		servers[i] = new ServerBasic(i, n);
	}

	int* jobsQueueSizes;
	jobsQueueSizes = new int[n]();

	// Dispatchers initialization
	DispatcherLocalServerState** dispatchers = new DispatcherLocalServerState * [m];

	// instantiate according to dispatcherType
	if (dispatcherType.compare("splittable_jsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableLocalServerStateJSQ(i, m, n, departure_parameters, seed);
		}
	}
	else
	{
		cout << "Dispatcher type received: " << dispatcherType << endl;
		assert("Unknown Dispatcher Type");
	}

	// Communication channels
	SizeTimestampPair* communicationChannels = new SizeTimestampPair[n]();
	for (int i = 0; i < n; ++i)
	{
		communicationChannels[i].first = 0;
		communicationChannels[i].second = NULL_TIMESTAMP;
	}

	// Main loop
	for (int t = 1; t < time + 1; ++t)
	{

		// collect queue sizes from servers
		for (int srv = 0; srv < n; ++srv)
		{
			jobsQueueSizes[srv] = servers[srv]->getQueueSize();
		}

		// logger: print status every 1% of time
		int curr_sum = 0; for (int i = 0; i < n; i++) curr_sum += jobsQueueSizes[i];
		logger.printStatus(t, time / 100, curr_sum);

		// logger: collect qls
		logger.updateMQL(jobsQueueSizes, n, t);

		// communication
		for (int disp = 0; disp < m; ++disp)
		{
			// clear communication channel
			for (int srv = 0; srv < n; ++srv)
			{
				communicationChannels[srv].first = 0;
				communicationChannels[srv].second = NULL_TIMESTAMP;
			}

			// randomly pick d out of n (no repeats)
			vector<int> sampled_queues;
			while (sampled_queues.size() < d)
			{
				int random_queue = unifDistribution(local_generator);
				if (find(sampled_queues.begin(), sampled_queues.end(), random_queue) == sampled_queues.end())
				{
					sampled_queues.push_back(random_queue);
				}
			}

			// update dispatcher local state
			for (int i = 0; i < d; ++i)
			{
				int srv = sampled_queues[i];
				servers[srv]->getInformation(communicationChannels, t);
				dispatchers[disp]->updateLocalState(communicationChannels);
			}
		}

		// dispatch arrivals
		int* tempJobsQueueSizes = new int[n];
		for (int disp = 0; disp < m; ++disp)
		{
			for (int srv = 0; srv < n; ++srv)
			{
				tempJobsQueueSizes[srv] = jobsQueueSizes[srv];
			}
			vector<DestinationJobsPairs> djPairs = dispatchers[disp]->getDestinations(rs.arrivals(disp));
			for (int j = 0; j < (int)djPairs.size(); ++j)
			{
				int destination = djPairs[j].first;
				int destinationJobs = djPairs[j].second;

				// send jobs to destination
				servers[destination]->insertJobs(destinationJobs, t);

				if (true) // syncWithDestinations
				{
					// add destinationJobs to destination at server local state
					tempJobsQueueSizes[destination] += destinationJobs;
					SizeTimestampPair pair = { tempJobsQueueSizes[destination] , t };
					dispatchers[disp]->updateLocalStateEntry(destination, pair);
				}
			}
		}
		delete[] tempJobsQueueSizes;

		// departures
		for (int i = 0; i < n; ++i)
		{
			int departures = rs.departures(i);
			vector<int> jobCreationTimes = servers[i]->serveJobsAndReturnCreationTimes(departures);
			for (int j = 0; j < (int)jobCreationTimes.size(); ++j)
			{
				// logger: collect jct
				logger.updateJCT(jobCreationTimes[j], t);
			}
		}

	}

	// logger: write results to a file
	logger.writeResultsToFiles(fn, time);

	// release memory

	delete[] jobsQueueSizes;

	for (int i = 0; i < n; i++)
	{
		delete servers[i];
	}
	delete[] servers;

	for (int i = 0; i < m; i++)
	{
		delete dispatchers[i];
	}
	delete[] dispatchers;

	delete[] communicationChannels;

}

void lsq_sample_het(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, int d, string dispatcherType)
{ // TODO

	// result file names
	string fn;
	fn.append(dispatcherType);
	fn.append("_lsq_");
	fn.append(to_string(d));
	fn.append("_chs_");
	fn.append(to_string(n));
	fn.append("_ser_");
	fn.append(to_string(m));
	fn.append("_dis_");
	fn.append(returnLoadString(arrival_parameters, m, departure_parameters, n));
	fn.append("_load_");
	fn.append(to_string(seed));
	fn.append("_seed_");

	// we use randomization here
	default_random_engine local_generator(seed);

	// for power of d - server randomization
	vector<double> auxilary_rinit_vec(departure_parameters, departure_parameters + n);
	discrete_distribution<int> server_sampling_distribution(auxilary_rinit_vec.begin(), auxilary_rinit_vec.end());
		
	// collect results
	Logger logger(time, fn_resolution);

	// arrivals and departures
	RandomProcesses rs(seed, m, arrival_parameters, n, departure_parameters);

	// servers initialization
	ServerBasic** servers = new ServerBasic * [n];
	for (int i = 0; i < n; i++)
	{
		servers[i] = new ServerBasic(i, n);
	}

	int* jobsQueueSizes;
	jobsQueueSizes = new int[n]();

	// Dispatchers initialization
	DispatcherLocalServerState** dispatchers = new DispatcherLocalServerState * [m];

	// instantiate according to dispatcherType
	if (dispatcherType.compare("splittable_hjfsq") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			dispatchers[i] = new DispatcherSplittableLocalServerStateHJFSQ(i, m, n, departure_parameters, seed);
		}
	}
	else
	{
		cout << "Dispatcher type received: " << dispatcherType << endl;
		assert("Unknown Dispatcher Type");
	}

	// Communication channels
	SizeTimestampPair* communicationChannels = new SizeTimestampPair[n]();
	for (int i = 0; i < n; ++i)
	{
		communicationChannels[i].first = 0;
		communicationChannels[i].second = NULL_TIMESTAMP;
	}

	// Main loop
	for (int t = 1; t < time + 1; ++t)
	{

		// collect queue sizes from servers
		for (int srv = 0; srv < n; ++srv)
		{
			jobsQueueSizes[srv] = servers[srv]->getQueueSize();
		}

		// logger: print status every 1% of time
		int curr_sum = 0; for (int i = 0; i < n; i++) curr_sum += jobsQueueSizes[i];
		logger.printStatus(t, time / 100, curr_sum);

		// logger: collect qls
		logger.updateMQL(jobsQueueSizes, n, t);

		// communication
		for (int disp = 0; disp < m; ++disp)
		{
			// clear communication channel
			for (int srv = 0; srv < n; ++srv)
			{
				communicationChannels[srv].first = 0;
				communicationChannels[srv].second = NULL_TIMESTAMP;
			}

			// weighted-randomly pick d out of n (no repeats) 
			vector<int> sampled_queues;
			while (sampled_queues.size() < d)
			{
				int random_queue = server_sampling_distribution(local_generator);
				if (find(sampled_queues.begin(), sampled_queues.end(), random_queue) == sampled_queues.end())
				{
					sampled_queues.push_back(random_queue);
				}
			}

			// update dispatcher local state
			for (int i = 0; i < d; ++i)
			{
				int srv = sampled_queues[i];
				servers[srv]->getInformation(communicationChannels, t);
				dispatchers[disp]->updateLocalState(communicationChannels);
			}
		}

		// dispatch arrivals
		int* tempJobsQueueSizes = new int[n];
		for (int disp = 0; disp < m; ++disp)
		{
			for (int srv = 0; srv < n; ++srv)
			{
				tempJobsQueueSizes[srv] = jobsQueueSizes[srv];
			}
			vector<DestinationJobsPairs> djPairs = dispatchers[disp]->getDestinations(rs.arrivals(disp));
			for (int j = 0; j < (int)djPairs.size(); ++j)
			{
				int destination = djPairs[j].first;
				int destinationJobs = djPairs[j].second;

				// send jobs to destination
				servers[destination]->insertJobs(destinationJobs, t);

				if (true) // syncWithDestinations
				{
					// add destinationJobs to destination at server local state
					tempJobsQueueSizes[destination] += destinationJobs;
					SizeTimestampPair pair = { tempJobsQueueSizes[destination] , t };
					dispatchers[disp]->updateLocalStateEntry(destination, pair);
				}
			}
		}
		delete[] tempJobsQueueSizes;

		// departures
		for (int i = 0; i < n; ++i)
		{
			int departures = rs.departures(i);
			vector<int> jobCreationTimes = servers[i]->serveJobsAndReturnCreationTimes(departures);
			for (int j = 0; j < (int)jobCreationTimes.size(); ++j)
			{
				// logger: collect jct
				logger.updateJCT(jobCreationTimes[j], t);
			}
		}

	}

	// logger: write results to a file
	logger.writeResultsToFiles(fn, time);

	// release memory

	delete[] jobsQueueSizes;

	for (int i = 0; i < n; i++)
	{
		delete servers[i];
	}
	delete[] servers;

	for (int i = 0; i < m; i++)
	{
		delete dispatchers[i];
	}
	delete[] dispatchers;

	delete[] communicationChannels;

}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////