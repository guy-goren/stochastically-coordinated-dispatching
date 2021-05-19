#include "RandomProcesses.h"

RandomProcesses::RandomProcesses(int seed, int m, double* arrival_parameters, int n, double* departure_parameters) : 
	n_(n), 
	m_(m), 
	generator_(seed)
{
	g_distribution_ = new geometric_distribution<int>*[n];
	for (int i = 0; i < n; ++i)
	{
		g_distribution_[i] = new geometric_distribution<int>(1.0 / (departure_parameters[i] + 1.0));
	}

	p_distribution_ = new poisson_distribution<int>*[m];
	for (int i = 0; i < m; ++i)
	{
		p_distribution_[i] = new poisson_distribution<int>(arrival_parameters[i]);
	}
}

RandomProcesses::~RandomProcesses()
{
	for (int i = 0; i < n_; ++i)
	{
		delete g_distribution_[i];
	}
	delete[] g_distribution_;

	for (int i = 0; i < m_; ++i)
	{
		delete p_distribution_[i];
	}
	delete[] p_distribution_;
}

int RandomProcesses::arrivals(int i)
{
	return (*p_distribution_[i])(generator_);
}

int RandomProcesses::departures(int i)
{
	return (*g_distribution_[i])(generator_);
}
