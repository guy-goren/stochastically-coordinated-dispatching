#include "defs.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class Dispatcher 
{

protected:

	int m_dispatcherID;
	int m_NumOfDispatchers;
	int m_NumberOfServers;

	double* m_ExpectedServiceCapacity;

	default_random_engine m_random_number_engine;

	// for fast WL computations
	int* m_wl_sorted_indices;
	double* m_wl_normalized_queue_lengths;

public:

	double ComputeWaterLevel(const int* queueLenghts, int EstimatedTotalArrivals);
	double ComputeWaterLevelFast(const int* queueLenghts, int EstimatedTotalArrivals);

	virtual vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int batchArrival)
	{ 
		return vector<DestinationJobsPairs>();
	}

	Dispatcher(int id, int m, int n, double* esc, int seed);
	virtual ~Dispatcher();

};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableFullServerState : public Dispatcher
{
private:

public:

	DispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerState();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class SlowHeterogenousDispatcherSplittableFullServerState : public Dispatcher
{
private:

public:

	SlowHeterogenousDispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed);
	~SlowHeterogenousDispatcherSplittableFullServerState();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class FastHeterogenousDispatcherSplittableFullServerState : public Dispatcher
{
private:

	// auxilary arrays
	int* m_sorted_indices;
	double* m_normalized_queue_lengths;
	double* m_probabilities;
	int* m_probabilities_indices;
	int* m_sent_wrt_probabilities;
	double* m_calculation_aux;

public:

	FastHeterogenousDispatcherSplittableFullServerState(int id, int m, int n, double* s, int seed);
	~FastHeterogenousDispatcherSplittableFullServerState();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableFullServerStateJSQ : public Dispatcher
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableFullServerStateJSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerStateJSQ();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableFullServerStateHJSQ : public Dispatcher
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableFullServerStateHJSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerStateHJSQ();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


class DispatcherSplittableFullServerStateJFSQ : public Dispatcher
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableFullServerStateJFSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerStateJFSQ();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableFullServerStateHJFSQ : public Dispatcher
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableFullServerStateHJFSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerStateHJFSQ();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableFullServerStateFastHJFSQ : public Dispatcher
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableFullServerStateFastHJFSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableFullServerStateFastHJFSQ();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherJIQ : public Dispatcher
{

public:

	DispatcherJIQ(int id, int m, int n, double* s, int seed);
	~DispatcherJIQ();

	virtual vector<DestinationJobsPairs> getDestinations(int arrivals) = 0;
	virtual void updateIdleQueue(int index) = 0;
};

////////////////////////////////////////////////////////////////////

class DispatcherSplittableJIQ : public DispatcherJIQ
{
private:

	uniform_int_distribution<int> m_unifDistribution;

	vector<int> m_idle_queues;

public:

	DispatcherSplittableJIQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableJIQ();

	vector<DestinationJobsPairs> getDestinations(int arrivals);
	void updateIdleQueue(int index);
};

////////////////////////////////////////////////////////////////////

class DispatcherSplittableHetJFIQ : public DispatcherJIQ
{
private:

	uniform_int_distribution<int> m_unifDistribution;
	discrete_distribution<int>* m_server_sampling_distribution;

	priority_queue<pair<pair<int, double>, int>, vector<pair<pair<int, double>, int>>, CustomCompareHJFIQ> m_idle_queues_pq;
	bool* m_idle_queues;

public:

	DispatcherSplittableHetJFIQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableHetJFIQ();

	vector<DestinationJobsPairs> getDestinations(int arrivals);
	void updateIdleQueue(int index);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherLocalServerState : public Dispatcher
{
protected:

	SizeTimestampPair* m_localState;

public:

	DispatcherLocalServerState(int id, int m, int n, double* s, int seed);
	virtual ~DispatcherLocalServerState();

	virtual vector<DestinationJobsPairs> getDestinations(int arrivals)
	{
		return vector<DestinationJobsPairs>();
	}

	void updateLocalStateEntry(int server, SizeTimestampPair pair);
	void updateLocalState(const SizeTimestampPair* sizeTimestampPairs);
	void getLocalState(SizeTimestampPair* sizeTimestampPairs);
};

////////////////////////////////////////////////////////////////////

class DispatcherSplittableLocalServerStateJSQ : public DispatcherLocalServerState
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableLocalServerStateJSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableLocalServerStateJSQ();

	vector<DestinationJobsPairs> getDestinations(int arrivals);
};

////////////////////////////////////////////////////////////////////

class DispatcherSplittableLocalServerStateHJFSQ : public DispatcherLocalServerState
{
private:

	uniform_int_distribution<int> m_unifDistribution;

public:

	DispatcherSplittableLocalServerStateHJFSQ(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableLocalServerStateHJFSQ();

	vector<DestinationJobsPairs> getDestinations(int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherTWF
{

protected:

	int m_dispatcherID;
	int m_NumOfDispatchers;
	int m_NumberOfServers;

	default_random_engine m_random_number_engine;

public:

	double ComputeWaterLevel(const int* queueLenghts, int EstimatedTotalArrivals);

	virtual vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int batchArrival)
	{
		return vector<DestinationJobsPairs>();
	}

	DispatcherTWF(int id, int m, int n, double* s, int seed);
	virtual ~DispatcherTWF();

};

////////////////////////////////////////////////////////////////////

class DispatcherTWFSplittableFullServerState : public DispatcherTWF
{
private:

public:

	DispatcherTWFSplittableFullServerState(int id, int m, int n, double* s, int seed);
	~DispatcherTWFSplittableFullServerState();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int batchArrival);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class DispatcherSplittableWeightedRandom : public Dispatcher
{
private:

	discrete_distribution<int>* m_server_sampling_distribution;

public:

	DispatcherSplittableWeightedRandom(int id, int m, int n, double* s, int seed);
	~DispatcherSplittableWeightedRandom();

	vector<DestinationJobsPairs> getDestinations(const int* queueLengths, int arrivals);
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
