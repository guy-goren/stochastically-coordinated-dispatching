#include "dispatcher.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

Dispatcher::Dispatcher(int id, int m, int n, double* esc, int seed) :
	m_dispatcherID(id),
	m_NumOfDispatchers(m),
	m_NumberOfServers(n),
	m_random_number_engine(seed + m*id)
{
	m_ExpectedServiceCapacity = new double[n];
	for (int i = 0; i < n; i++)
	{
		m_ExpectedServiceCapacity[i] = esc[i];
	}

	m_wl_sorted_indices = new int[m_NumberOfServers];
	m_wl_normalized_queue_lengths = new double[m_NumberOfServers];
}

Dispatcher::~Dispatcher()
{
	delete[] m_ExpectedServiceCapacity;

	delete[] m_wl_sorted_indices;
	delete[] m_wl_normalized_queue_lengths;
}

double Dispatcher::ComputeWaterLevel(const int * queueLenghts, int EstimatedTotalArrivals)
{
	
	// keep normalized height (i.e., q_s/mu_s) in a "mean heap" priority queue
	priority_queue<pair<int, double>, vector<pair<int,double>>, CustomCompareWL> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push( make_pair(queueLenghts[i] , m_ExpectedServiceCapacity[i]) );
	}

	// init
	double mu_tot = 0;
	double remainning_jobs = (double)EstimatedTotalArrivals;
	pair<int, double> curr_q = pq.top();
	double waterLevel = curr_q.first / curr_q.second;
	double heigt_to_fill = 0;
	
	// "pour" jobs to queues 
	while (remainning_jobs > 0) {
		if (pq.empty()) {
			return waterLevel + remainning_jobs / mu_tot;
		}
		else {
			mu_tot += curr_q.second;
			pq.pop();
			if (pq.empty()) {
				return  waterLevel + remainning_jobs / mu_tot;
			}
			curr_q = pq.top();
			heigt_to_fill = (curr_q.first / curr_q.second) - waterLevel;

			if (heigt_to_fill * mu_tot < remainning_jobs) {
				waterLevel += heigt_to_fill;
				remainning_jobs -= heigt_to_fill * mu_tot;
			}
			else {
				return waterLevel + remainning_jobs / mu_tot;
			}
		}
	}

	// should never reach here
	assert("Error in Dispatcher::ComputeWaterLevel: reached end of function");

	return -1;
}

double Dispatcher::ComputeWaterLevelFast(const int* queueLenghts, int EstimatedTotalArrivals)
{

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		m_wl_sorted_indices[i] = i;
		m_wl_normalized_queue_lengths[i] = queueLenghts[i] / m_ExpectedServiceCapacity[i];
	}

	// after this, sorted_indices contains the indices of sorted values in "normalized_queue_lengths"
	sort(m_wl_sorted_indices, m_wl_sorted_indices + m_NumberOfServers, sort_indices(m_wl_normalized_queue_lengths));

	double remaining_jobs = (double)EstimatedTotalArrivals;

	double total_witdh = 0;
	double water_level = m_wl_normalized_queue_lengths[m_wl_sorted_indices[0]];

	for (int i = 0; i < m_NumberOfServers; ++i)
	{

		total_witdh += m_ExpectedServiceCapacity[m_wl_sorted_indices[i]];

		if (i == m_NumberOfServers - 1)
		{
			return water_level + remaining_jobs / total_witdh;
		}

		double t = m_wl_normalized_queue_lengths[m_wl_sorted_indices[i + 1]];
		double delta = t - water_level;
		double to_dispatch = delta * total_witdh;
		if (to_dispatch >= remaining_jobs)
		{
			return water_level + remaining_jobs / total_witdh;
		}

		remaining_jobs -= to_dispatch;
		water_level += delta;

	}

	// should never reach here
	assert("Error in Dispatcher::ComputeWaterLevel: reached end of function");

	return -1;
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerState::DispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
}

DispatcherSplittableFullServerState::~DispatcherSplittableFullServerState()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerState::getDestinations(const int * queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{		
		return djPairs;
	}

	// estimate total arrivals (at all dispatchers) and compute water-level
	// beyond this point, it must hold that "EstimatedTotalArrivals > 1"
	int EstimatedTotalArrivals = arrivals * m_NumOfDispatchers;
	double waterLevel = ComputeWaterLevel(queueLengths, EstimatedTotalArrivals);

	// sort according to (2q_s+1) / mu_s
	priority_queue<pair<int, double>, vector<pair<int, double>>, CustomCompareGD> pq;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push( make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]) );
	}

	// init
	vector<pair<int, double>> curr_set;	
	double opt_Lambda0(0);
	double Lambda0_nom = -2.0*(EstimatedTotalArrivals - 1.0);
	double Lambda0_denom(0);
	double val1(0);
	double val2(0);
	double curr_val;

	pair<int, double> curr_q;

	// initilize comparison values that are "infty" in the psuedocode
	curr_q = pq.top();
	double min_serv_val_for_feasibilty_test = 2 * waterLevel - (2.0 * curr_q.first + 1) / curr_q.second;
	double Lambda0 = (-2.0 * (EstimatedTotalArrivals - 1.0) + 2 * (curr_q.second * waterLevel - curr_q.first) - 1) / curr_q.second;
	double opt_val = curr_q.second / (4.0 * (EstimatedTotalArrivals - 1.0)) * pow( Lambda0, 2) - 
		pow(2 * (curr_q.first - curr_q.second * waterLevel) + 1, 2) / (4 * curr_q.second * (EstimatedTotalArrivals - 1.0));

	while (!pq.empty() ) {

		curr_q = pq.top();
		curr_set.push_back(curr_q);
		pq.pop();
		
		// compute current Lambda_0 - for the current set of servers
		Lambda0_nom += 2 * (curr_q.second * waterLevel - curr_q.first) - 1;
		Lambda0_denom += curr_q.second;
		Lambda0 = Lambda0_nom / Lambda0_denom;

		// compare only worst case - the server with the lowest value
		min_serv_val_for_feasibilty_test = std::min(min_serv_val_for_feasibilty_test, 2 * waterLevel - (2.0 * curr_q.first + 1) / curr_q.second);
		
		if (min_serv_val_for_feasibilty_test - Lambda0 < 0) { continue; } // a probability<0 exists, so this solution is infeasible
		val1 += curr_q.second / (4.0 * (EstimatedTotalArrivals - 1.0));
		val2 += pow(2 * (curr_q.first - curr_q.second * waterLevel) + 1, 2) / (4 * curr_q.second * (EstimatedTotalArrivals - 1.0));
		curr_val = val1*pow(Lambda0,2) - val2;

		if (curr_val < opt_val) { // is this the best solution so far
			opt_val = curr_val;
			opt_Lambda0 = Lambda0;
		}

	}

	//compute the optimal send probabilities based on opt_Lambda
	vector<double> probabilities(m_NumberOfServers, 0);

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		// compute its probability and remember it
		double temp = (-2 * (queueLengths[i] - m_ExpectedServiceCapacity[i] * waterLevel) - 1 - m_ExpectedServiceCapacity[i] * opt_Lambda0) / (2.0 * (EstimatedTotalArrivals - 1.0));
		double probability = std::max(0.0, temp);
		probabilities[i] = probability;
	}

	// create the discrete distribution
	discrete_distribution<int> destination_distribution(probabilities.begin(), probabilities.end());

	// randomize the destinations according to probabilities
	for (int i = 0; i < arrivals; ++i)
	{
		// get the random destination
		int destination = destination_distribution(m_random_number_engine);	
		
		// create and admit the result
		DestinationJobsPairs djPair(destination, 1);
		djPairs.push_back(djPair);
	}

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SlowHeterogenousDispatcherSplittableFullServerState::SlowHeterogenousDispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
}

SlowHeterogenousDispatcherSplittableFullServerState::~SlowHeterogenousDispatcherSplittableFullServerState()
{
}

vector<DestinationJobsPairs> SlowHeterogenousDispatcherSplittableFullServerState::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// estimate total arrivals (at all dispatchers) and compute water-level
	// beyond this point, it must hold that "EstimatedTotalArrivals > 1"
	int EstimatedTotalArrivals = arrivals * m_NumOfDispatchers;
	double waterLevel = ComputeWaterLevel(queueLengths, EstimatedTotalArrivals);

	// sort according to (2q_s+1) / mu_s
	priority_queue<pair<int, double>, vector<pair<int, double>>, CustomCompareGD> pq;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]));
	}

	// init
	vector<pair<int, double>> curr_set;
	vector<double> curr_probabilities;
	for (int i = 0; i < m_NumberOfServers; i++)
	{
		curr_probabilities.push_back(0.0);
	}
	double Lambda0(0);
	double Lambda0_nom = (0);
	double Lambda0_denom(0);
	double opt_Lambda0(0);
	double probability(0);
	double curr_val;
	double opt_val(0);

	pair<int, double> curr_q;

	// initilize comparison values that are "infty" in the psuedocode
	curr_q = pq.top();
	Lambda0 = (-2.0 * (EstimatedTotalArrivals - 1.0) + 2 * (curr_q.second * waterLevel - curr_q.first) - 1) / curr_q.second;
	opt_val = curr_q.second / (4.0 * (EstimatedTotalArrivals - 1.0)) * pow(Lambda0, 2) -
		pow(2 * (curr_q.first - curr_q.second * waterLevel) + 1, 2) / (4 * curr_q.second * (EstimatedTotalArrivals - 1.0));

	while (!pq.empty()) {

		curr_q = pq.top();
		curr_set.push_back(curr_q);
		pq.pop();

		// compute current Lambda_0 - for the current set of servers
		Lambda0_nom = -2.0 * (EstimatedTotalArrivals - 1.0);
		Lambda0_denom = (0);
		for (int i = 0; i < curr_set.size(); i++)
		{
			Lambda0_nom += 2 * (curr_set.at(i).second * waterLevel - curr_set.at(i).first) - 1;
			Lambda0_denom += curr_set.at(i).second;
		}
		Lambda0 = (Lambda0_denom > 0) ? Lambda0_nom / Lambda0_denom : 0;

		// check feasibility of the current probabilities (p>0)...
		//    ...and compute the value of the objective function for the current set
		curr_val = 0.0;
		for (int i = 0; i < curr_set.size(); i++)
		{
			probability = (-2 * (curr_set.at(i).first - curr_set.at(i).second * waterLevel) - 1 - curr_set.at(i).second * Lambda0 )
				/ (2.0 * (EstimatedTotalArrivals - 1.0));
			if (probability < 0) // infeasible probability
			{
				goto OUTER_CONTINUE;
			}
			curr_val += ((EstimatedTotalArrivals - 1.0) * pow(probability, 2) + (2 * (curr_set.at(i).first - curr_set.at(i).second * waterLevel) + 1) * probability) / curr_set.at(i).second;
		}
		OUTER_CONTINUE: continue; // continue to the next set due to an infeasible probability of the current
				
		if (curr_val < opt_val) { // is this the best solution so far
			opt_val = curr_val;
			opt_Lambda0 = Lambda0;
		}
	}

	//compute the optimal send probabilities based on opt_Lambda
	vector<double> probabilities(m_NumberOfServers, 0);

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		// compute its probability and remember it
		double temp = (-2 * (queueLengths[i] - m_ExpectedServiceCapacity[i] * waterLevel) - 1 - m_ExpectedServiceCapacity[i] * opt_Lambda0) / (2.0 * (EstimatedTotalArrivals - 1.0));
		probability = std::max(0.0, temp);
		probabilities[i] = probability;
	}

	// create the discrete distribution
	discrete_distribution<int> destination_distribution(probabilities.begin(), probabilities.end());

	// randomize the destinations according to probabilities
	for (int i = 0; i < arrivals; ++i)
	{
		// get the random destination
		int destination = destination_distribution(m_random_number_engine);

		// create and admit the result
		DestinationJobsPairs djPair(destination, 1);
		djPairs.push_back(djPair);
	}

	return djPairs;
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

FastHeterogenousDispatcherSplittableFullServerState::FastHeterogenousDispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
	// prepare auxilary arrays
	m_sorted_indices = new int[m_NumberOfServers];
	m_normalized_queue_lengths = new double[m_NumberOfServers];
	m_probabilities = new double[m_NumberOfServers]();
	m_sent_wrt_probabilities = new int[m_NumberOfServers];
	m_calculation_aux = new double[m_NumberOfServers];
	m_probabilities_indices = new int[m_NumberOfServers];
}

FastHeterogenousDispatcherSplittableFullServerState::~FastHeterogenousDispatcherSplittableFullServerState()
{
	delete[] m_sorted_indices;
	delete[] m_normalized_queue_lengths;
	delete[] m_probabilities;
	delete[] m_sent_wrt_probabilities;
	delete[] m_calculation_aux;
	delete[] m_probabilities_indices;
}

vector<DestinationJobsPairs> FastHeterogenousDispatcherSplittableFullServerState::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// estimate total arrivals (at all dispatchers) and compute water-level
	// beyond this point, it must hold that "EstimatedTotalArrivals > 1"
	int EstimatedTotalArrivals = arrivals * m_NumOfDispatchers;
	double waterLevel = ComputeWaterLevelFast(queueLengths, EstimatedTotalArrivals);

	// initialize auxilary arrays
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		m_sorted_indices[i] = i;
		m_normalized_queue_lengths[i] = (2.0 * queueLengths[i] + 1.0) / m_ExpectedServiceCapacity[i];
	}

	// after this, sorted_indices contains the indices of sorted values in "normalized_queue_lengths"
	sort(m_sorted_indices, m_sorted_indices + m_NumberOfServers, sort_indices(m_normalized_queue_lengths));

	double v_1 = 0;
	double v_2 = 0;

	double lambda_0_n = -2.0 * (EstimatedTotalArrivals - 1.0);
	double lambda_0_d = 0;

	double x = 2 * waterLevel;
	double best_val = LONG_MAX;
	double best_lambda_0 = LONG_MAX;

	double EstimatedTotalArrivalsMinusOne = EstimatedTotalArrivals - 1.0;

	// m_calculation_aux
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		m_calculation_aux[i] = 2 * (queueLengths[i] - m_ExpectedServiceCapacity[i] * waterLevel) + 1;
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{

		int current_index = m_sorted_indices[i];

		double t = m_normalized_queue_lengths[current_index];

		lambda_0_n -= m_calculation_aux[current_index];
		lambda_0_d += m_ExpectedServiceCapacity[current_index];

		double lambda_0 = lambda_0_n / lambda_0_d;

		double candidate_x = 2 * waterLevel - m_normalized_queue_lengths[current_index];
		x = (x < candidate_x) ? x : candidate_x;

		if (x - lambda_0 < 0)
		{
			continue;
		}

		v_1 += m_ExpectedServiceCapacity[current_index] / (4 * EstimatedTotalArrivalsMinusOne);

		double auxilary_v2_1 = m_calculation_aux[current_index];
		double auxilary_v2_2 = 4 * m_ExpectedServiceCapacity[current_index] * EstimatedTotalArrivalsMinusOne;
		v_2 += auxilary_v2_1 * auxilary_v2_1 / auxilary_v2_2;

		double val = v_1 * lambda_0 * lambda_0 - v_2;
		if (val < best_val)
		{
			best_val = val;
			best_lambda_0 = lambda_0;
		}

	}

	// compute the probabilities
	int num_positive_probabilities = 0;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		double temp = (-m_calculation_aux[i] - m_ExpectedServiceCapacity[i] * best_lambda_0) / (2 * EstimatedTotalArrivalsMinusOne);
		if (temp > 0)
		{
			m_probabilities[num_positive_probabilities] = temp;
			m_probabilities_indices[num_positive_probabilities] = i;
			num_positive_probabilities += 1;
		}
	}

	// create the discrete distribution
	discrete_distribution<int> destination_distribution(m_probabilities, m_probabilities + num_positive_probabilities);

	// randomize the destinations according to probabilities
	memset(m_sent_wrt_probabilities, 0, sizeof(int) * m_NumberOfServers);
	for (int i = 0; i < arrivals; ++i)
	{
		// get the random destination
		int destination = destination_distribution(m_random_number_engine);
		m_sent_wrt_probabilities[m_probabilities_indices[destination]] += 1;
	}
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (m_sent_wrt_probabilities[i] > 0)
		{
			// create and admit the result
			DestinationJobsPairs djPair(i, m_sent_wrt_probabilities[i]);
			djPairs.push_back(djPair);
		}
	}

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerStateJSQ::DispatcherSplittableFullServerStateJSQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableFullServerStateJSQ::~DispatcherSplittableFullServerStateJSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerStateJSQ::getDestinations(const int * queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// compute the "water level"
	priority_queue<int, vector<int>, greater<int>> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(queueLengths[i]);
	}

	for (int i = 0; i < arrivals; ++i)
	{
		pq.push(pq.top() + 1);
		pq.pop();
	}

	int water_level = pq.top();

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] < water_level)
		{
			int sent = (water_level - queueLengths[i]);			
			DestinationJobsPairs djPair(i, sent);
			djPairs.push_back(djPair);
			cumulativeSent += sent;
		}
	}

	// no remainder
	if (cumulativeSent == arrivals)
	{
		return djPairs;
	}

	// send the remainder randomly to shortest queues
	vector<int> leqWaterLevel;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] <= water_level)
		{
			leqWaterLevel.push_back(i);
		}
	}

	for (int i = 0; i < arrivals - cumulativeSent; ++i)
	{
		// chose a queue
		int randomIndex = m_unifDistribution(m_random_number_engine) % leqWaterLevel.size();
		int destination = leqWaterLevel[randomIndex];

		// send a single job to that queue
		DestinationJobsPairs djPair(destination, 1);
		djPairs.push_back(djPair);

		// erase it from candidates
		leqWaterLevel.erase(leqWaterLevel.begin() + randomIndex);
	}

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerStateHJSQ::DispatcherSplittableFullServerStateHJSQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableFullServerStateHJSQ::~DispatcherSplittableFullServerStateHJSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerStateHJSQ::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// sort according to queue length and power
	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareHJSQ> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{		
		pq.push(make_pair(make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]), i));
	}

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	int* send = new int[m_NumberOfServers]();

	while (cumulativeSent < arrivals)
	{
		int qle = pq.top().first.first;
		double esr = pq.top().first.second;
		int dst = pq.top().second;

		send[dst] += 1;
		cumulativeSent += 1;
		
		pq.pop();
		pq.push(make_pair(make_pair(qle + 1, esr), dst));
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;
	
	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerStateJFSQ::DispatcherSplittableFullServerStateJFSQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableFullServerStateJFSQ::~DispatcherSplittableFullServerStateJFSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerStateJFSQ::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// sort according to queue length and power
	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareJFSQ> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(make_pair(make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]), i));
	}

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	int* send = new int[m_NumberOfServers]();

	while (cumulativeSent < arrivals)
	{
		int qle = pq.top().first.first;
		double esr = pq.top().first.second;
		int dst = pq.top().second;

		send[dst] += 1;
		cumulativeSent += 1;

		pq.pop();
		pq.push(make_pair(make_pair(qle + 1, esr), dst));
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerStateHJFSQ::DispatcherSplittableFullServerStateHJFSQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableFullServerStateHJFSQ::~DispatcherSplittableFullServerStateHJFSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerStateHJFSQ::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// sort according to queue length and power
	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareHJFSQ> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(make_pair(make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]), i));
	}

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	int* send = new int[m_NumberOfServers]();

	while (cumulativeSent < arrivals)
	{
		int qle = pq.top().first.first;
		double esr = pq.top().first.second;
		int dst = pq.top().second;

		send[dst] += 1;
		cumulativeSent += 1;

		pq.pop();
		pq.push(make_pair(make_pair(qle + 1, esr), dst));
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableFullServerStateFastHJFSQ::DispatcherSplittableFullServerStateFastHJFSQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableFullServerStateFastHJFSQ::~DispatcherSplittableFullServerStateFastHJFSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableFullServerStateFastHJFSQ::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}


	if (arrivals == 1)
	{

		double min_val = queueLengths[0] / m_ExpectedServiceCapacity[0];
		int min_index = 0;

		for (int i = 1; i < m_NumberOfServers; ++i)
		{
			double val = queueLengths[i] / m_ExpectedServiceCapacity[i];
			if ((val < min_val) || ((val == min_val) && (m_ExpectedServiceCapacity[min_index] < m_ExpectedServiceCapacity[i])))
			{
				min_val = val;
				min_index = i;
			}
		}

		DestinationJobsPairs djPair(min_index, 1);
		djPairs.push_back(djPair);

		return djPairs;
	}


	double waterLevel = ComputeWaterLevelFast(queueLengths, arrivals);

	// sort according to queue length and power
	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareHJFSQ> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] / m_ExpectedServiceCapacity[i] <= waterLevel)
		{
			pq.push(make_pair(make_pair(queueLengths[i], m_ExpectedServiceCapacity[i]), i));
		}
	}

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	int* send = new int[m_NumberOfServers]();

	while (cumulativeSent < arrivals)
	{
		int qle = pq.top().first.first;
		double esr = pq.top().first.second;
		int dst = pq.top().second;

		send[dst] += 1;
		cumulativeSent += 1;

		pq.pop();
		pq.push(make_pair(make_pair(qle + 1, esr), dst));
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherJIQ::DispatcherJIQ(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
}

DispatcherJIQ::~DispatcherJIQ()
{
}

////////////////////////////////////////////////////////////////////

DispatcherSplittableJIQ::DispatcherSplittableJIQ(int id, int m, int n, double* s, int seed) :
	DispatcherJIQ(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableJIQ::~DispatcherSplittableJIQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableJIQ::getDestinations(int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// there are idle queues
	if (m_idle_queues.size() > 0)
	{
		int cumulativeSent = 0;
		int numIdles = (int)m_idle_queues.size();
		int idle_water_level = arrivals / numIdles;

		for (int i = 0; i < numIdles; ++i)
		{
			DestinationJobsPairs djPair(m_idle_queues[i], idle_water_level);
			djPairs.push_back(djPair);
			cumulativeSent += idle_water_level;
		}

		vector<int> idle_queues_remainder;
		while ((int)idle_queues_remainder.size() < arrivals - cumulativeSent)
		{
			int random_index = m_unifDistribution(m_random_number_engine) % m_idle_queues.size();
			idle_queues_remainder.push_back(m_idle_queues[random_index]);
			m_idle_queues.erase(m_idle_queues.begin() + random_index);
		}

		for (int i = 0; i < idle_queues_remainder.size(); ++i)
		{
			DestinationJobsPairs djPair(idle_queues_remainder[i], 1);
			djPairs.push_back(djPair);
			cumulativeSent += 1;
		}
	}
	else
	{
		for (int i = 0; i < arrivals; ++i)
		{
			int random_server = m_unifDistribution(m_random_number_engine) % m_NumberOfServers;
			DestinationJobsPairs djPair(random_server, 1);
			djPairs.push_back(djPair);
		}
	}

	m_idle_queues.clear();

	return djPairs;

}

void DispatcherSplittableJIQ::updateIdleQueue(int index)
{
	if (find(m_idle_queues.begin(), m_idle_queues.end(), index) == m_idle_queues.end())
	{
		m_idle_queues.push_back(index);
	}
}

////////////////////////////////////////////////////////////////////

DispatcherSplittableHetJFIQ::DispatcherSplittableHetJFIQ(int id, int m, int n, double* s, int seed) :
	DispatcherJIQ(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
	vector<double> vec(s, s + n);
	m_server_sampling_distribution = new discrete_distribution<int>(vec.begin(), vec.end());

	m_idle_queues = new bool[n]();
}

DispatcherSplittableHetJFIQ::~DispatcherSplittableHetJFIQ()
{
	delete m_server_sampling_distribution;
	delete[] m_idle_queues;
}

vector<DestinationJobsPairs> DispatcherSplittableHetJFIQ::getDestinations(int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	int cumulativeSent = 0;

	// there are idle queues
	bool has_idle_queues = false;

	while (m_idle_queues_pq.size() > 0)
	{
		m_idle_queues_pq.pop();
	}

	for (int i = 0; i < m_NumberOfServers; i++)
	{
		if (m_idle_queues[i])
		{
			has_idle_queues = true;
			m_idle_queues_pq.push(make_pair(make_pair(0, m_ExpectedServiceCapacity[i]), i));
		}
	}

	if (has_idle_queues)
	{

		int cumulativeSent = 0;

		int* send = new int[m_NumberOfServers]();

		while (cumulativeSent < arrivals)
		{
			int qle = m_idle_queues_pq.top().first.first;
			double esr = m_idle_queues_pq.top().first.second;
			int dst = m_idle_queues_pq.top().second;

			send[dst] += 1;
			cumulativeSent += 1;

			m_idle_queues_pq.pop();
			m_idle_queues_pq.push(make_pair(make_pair(qle + 1, esr), dst));

			m_idle_queues[dst] = false;
		}

		for (int i = 0; i < m_NumberOfServers; ++i)
		{
			if (send[i] > 0)
			{
				DestinationJobsPairs djPair(i, send[i]);
				djPairs.push_back(djPair);
			}
		}

		delete[] send;

	}
	else
	{ 
		for (int i = 0; i < arrivals; ++i)
		{
			int random_server = (*m_server_sampling_distribution)(m_random_number_engine);
			DestinationJobsPairs djPair(random_server, 1);
			djPairs.push_back(djPair);
		}
	}

	return djPairs;

}

void DispatcherSplittableHetJFIQ::updateIdleQueue(int index)
{
	m_idle_queues[index] = true;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherLocalServerState::DispatcherLocalServerState(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
	m_localState = new SizeTimestampPair[n]();
	for (int i = 0; i < n; ++i)
	{
		m_localState[i].first = 0;
		m_localState[i].second = NULL_TIMESTAMP;
	}
}

DispatcherLocalServerState::~DispatcherLocalServerState()
{
	delete[] m_localState;
}

void DispatcherLocalServerState::updateLocalStateEntry(int server, SizeTimestampPair pair)
{
	m_localState[server] = pair;
}

void DispatcherLocalServerState::updateLocalState(const SizeTimestampPair* sizeTimestampPairs)
{
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		// if the local timestamp is older update both queue length and its timestamp
		if (m_localState[i].second < sizeTimestampPairs[i].second)
		{
			m_localState[i] = sizeTimestampPairs[i];
		}
	}
}

void DispatcherLocalServerState::getLocalState(SizeTimestampPair* sizeTimestampPairs)
{
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		sizeTimestampPairs[i] = m_localState[i];
	}
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableLocalServerStateJSQ::DispatcherSplittableLocalServerStateJSQ(int id, int m, int n, double* s, int seed) :
	DispatcherLocalServerState(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableLocalServerStateJSQ::~DispatcherSplittableLocalServerStateJSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableLocalServerStateJSQ::getDestinations(int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// retreive local queue lengths from local state
	int* queueLengths = new int[m_NumberOfServers];
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		queueLengths[i] = m_localState[i].first;
	}

	// compute the "water level"
	priority_queue<int, vector<int>, greater<int>> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(queueLengths[i]);
	}

	for (int i = 0; i < arrivals; ++i)
	{
		pq.push(pq.top() + 1);
		pq.pop();
	}

	int water_level = pq.top();

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] < water_level)
		{
			int sent = (water_level - queueLengths[i]);
			DestinationJobsPairs djPair(i, sent);
			djPairs.push_back(djPair);
			cumulativeSent += sent;
		}
	}

	// no remainder
	if (cumulativeSent == arrivals)
	{
		delete[] queueLengths;
		return djPairs;
	}

	// send the remainder randomly to shortest queues
	vector<int> leqWaterLevel;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] <= water_level)
		{
			leqWaterLevel.push_back(i);
		}
	}

	for (int i = 0; i < arrivals - cumulativeSent; ++i)
	{
		// chose a queue
		int randomIndex = m_unifDistribution(m_random_number_engine) % leqWaterLevel.size();
		int destination = leqWaterLevel[randomIndex];

		// send a single job to that queue
		DestinationJobsPairs djPair(destination, 1);
		djPairs.push_back(djPair);

		// erase it from candidates
		leqWaterLevel.erase(leqWaterLevel.begin() + randomIndex);
	}

	delete[] queueLengths;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableLocalServerStateHJFSQ::DispatcherSplittableLocalServerStateHJFSQ(int id, int m, int n, double* s, int seed) :
	DispatcherLocalServerState(id, m, n, s, seed),
	m_unifDistribution(0, 100000000)
{
}

DispatcherSplittableLocalServerStateHJFSQ::~DispatcherSplittableLocalServerStateHJFSQ()
{
}

vector<DestinationJobsPairs> DispatcherSplittableLocalServerStateHJFSQ::getDestinations(int arrivals)
{ 

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// sort according to queue length and power
	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareHJFSQ> pq;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(make_pair(make_pair(m_localState[i].first, m_ExpectedServiceCapacity[i]), i));
	}

	// send jobs according to the "integer water level"
	int cumulativeSent = 0;

	int* send = new int[m_NumberOfServers]();

	while (cumulativeSent < arrivals)
	{
		int qle = pq.top().first.first;
		double esr = pq.top().first.second;
		int dst = pq.top().second;

		send[dst] += 1;
		cumulativeSent += 1;

		pq.pop();
		pq.push(make_pair(make_pair(qle + 1, esr), dst));
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherTWF::DispatcherTWF(int id, int m, int n, double* s, int seed) :
	m_dispatcherID(id),
	m_NumOfDispatchers(m),
	m_NumberOfServers(n),
	m_random_number_engine(seed + m * id)
{
}

DispatcherTWF::~DispatcherTWF()
{
}

double DispatcherTWF::ComputeWaterLevel(const int* queueLenghts, int EstimatedTotalArrivals)
{

	priority_queue<int, vector<int>, greater<int>> pq;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		pq.push(queueLenghts[i]);
	}

	int curr_wl_height = pq.top();
	pq.pop();
	int num_of_shortest_queues = 1;

	while ((!pq.empty()) && (curr_wl_height == pq.top()))
	{
		num_of_shortest_queues++;
		pq.pop();
	}

	int cumulative_sent = 0;

	while ((EstimatedTotalArrivals - cumulative_sent) / num_of_shortest_queues >= 1)
	{

		int heigt_to_fill = 1;

		if (!pq.empty())
		{
			heigt_to_fill = pq.top() - curr_wl_height;
		}

		heigt_to_fill = min(heigt_to_fill, ((EstimatedTotalArrivals - cumulative_sent) / num_of_shortest_queues));

		cumulative_sent += heigt_to_fill * num_of_shortest_queues;
		curr_wl_height += heigt_to_fill;

		while ((!pq.empty()) && (pq.top() == curr_wl_height))
		{
			pq.pop();
			num_of_shortest_queues++;
		}

	}

	double waterLevel = curr_wl_height + ((double)EstimatedTotalArrivals - (double)cumulative_sent) / (double)num_of_shortest_queues;

	return waterLevel;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherTWFSplittableFullServerState::DispatcherTWFSplittableFullServerState(int id, int m, int n, double* s, int seed) :
	DispatcherTWF(id, m, n, s, seed)
{
}

DispatcherTWFSplittableFullServerState::~DispatcherTWFSplittableFullServerState()
{
}

vector<DestinationJobsPairs> DispatcherTWFSplittableFullServerState::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	// estimate total arrivals (at all dispatchers) and compute water-level
	// beyond this point, it must hold that "EstimatedTotalArrivals > 1"
	int EstimatedTotalArrivals = arrivals * m_NumOfDispatchers;
	double waterLevel = ComputeWaterLevel(queueLengths, EstimatedTotalArrivals);

	// remember potential destinations and their amount
	vector<int> potential_destinations;
	int numServersUnderWL = 0;
	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (queueLengths[i] < waterLevel)
		{
			++numServersUnderWL;
			potential_destinations.push_back(i);
		}
	}

	// compute send probabilities to all candidates
	vector<double> probabilities;
	for (int i = 0; i < (int)potential_destinations.size(); ++i)
	{
		// computer weight of potential destination
		int potential_destination = potential_destinations[i];
		double belowWaterLevel = waterLevel - queueLengths[potential_destination];

		// compute its probability and remember it
		double probability = (belowWaterLevel - 1.0 / (double)numServersUnderWL) / ((double)EstimatedTotalArrivals - 1.0);

		// fix numerical error
		probability = probability > 0 ? probability : 0;

		probabilities.push_back(probability);
	}

	// create the discrete distribution
	discrete_distribution<int> destination_distribution(probabilities.begin(), probabilities.end());

	// randomize the destinations according to probabilities
	for (int i = 0; i < arrivals; ++i)
	{
		// get the random destination
		int destination_index = destination_distribution(m_random_number_engine);
		int destination = potential_destinations[destination_index];

		// create and admit the result
		DestinationJobsPairs djPair(destination, 1);
		djPairs.push_back(djPair);
	}

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DispatcherSplittableWeightedRandom::DispatcherSplittableWeightedRandom(int id, int m, int n, double* s, int seed) :
	Dispatcher(id, m, n, s, seed)
{
	vector<double> vec(s, s + n);
	m_server_sampling_distribution = new discrete_distribution<int>(vec.begin(), vec.end());
}

DispatcherSplittableWeightedRandom::~DispatcherSplittableWeightedRandom()
{
	delete m_server_sampling_distribution;
}

vector<DestinationJobsPairs> DispatcherSplittableWeightedRandom::getDestinations(const int* queueLengths, int arrivals)
{

	vector<DestinationJobsPairs> djPairs;

	// no arrivals? return empty vector
	if (arrivals == 0)
	{
		return djPairs;
	}

	int* send = new int[m_NumberOfServers]();

	int cumulativeSent = 0;
	while (cumulativeSent < arrivals)
	{
		int dst = (*m_server_sampling_distribution)(m_random_number_engine);

		send[dst] += 1;
		cumulativeSent += 1;
	}

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		if (send[i] > 0)
		{
			DestinationJobsPairs djPair(i, send[i]);
			djPairs.push_back(djPair);
		}
	}

	delete[] send;

	return djPairs;

}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
