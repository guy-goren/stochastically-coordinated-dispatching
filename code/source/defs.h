#ifndef DEFS_H_
#define DEFS_H_

#include <stdlib.h>
#include <iostream>
#include <random>
#include <queue>
#include <algorithm>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string>
#include <unordered_map>
#include <queue>
#include <numeric>
#include <algorithm>
#include <chrono>

using std::cout;
using std::cin;
using std::endl;

using std::find;
using std::abs;
using std::to_string;
using std::min_element;
using std::make_pair;
using std::accumulate;
using std::max_element;
using std::min;
using std::sort;

using std::unordered_map;
using std::queue;
using std::priority_queue;
using std::pair;
using std::greater;
using std::vector;
using std::string;
using std::ofstream;

using std::default_random_engine;

using std::poisson_distribution;
using std::geometric_distribution;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::discrete_distribution; 
using std::bernoulli_distribution;

struct job_t 
{
	int creation_time;
};

typedef pair<int, int> q_Qlen;
typedef pair<int, int> StoppingPoint_Index;
typedef pair<int, int> DestinationJobsPairs;
typedef pair<int, int> SizeTimestampPair;

#define NULL_TIMESTAMP -1

struct CustomCompareWL
{
	bool operator()(const pair<int, double> lhs, const pair<int, double> rhs)
	{
		return (lhs.first / lhs.second) > (rhs.first / rhs.second);
	}
};

struct CustomCompareGD
{
	bool operator()(const pair<int, double> lhs, const pair<int, double> rhs)
	{
		return ((2.0 * lhs.first + 1) / lhs.second) > ((2.0 * rhs.first + 1) / rhs.second);
	}
};

struct CustomCompareHJSQ
{
	bool operator()(const pair<pair<int, double>, int> lhs, const pair<pair<int, double>, int> rhs)
	{
		if ((lhs.first.first / lhs.first.second) != (rhs.first.first / rhs.first.second))
		{
			return (lhs.first.first / lhs.first.second) > (rhs.first.first / rhs.first.second);
		}
		else
		{
			return (bool)(rand() % 2);
		}
	}
};

struct CustomCompareJFSQ
{
	bool operator()(const pair<pair<int, double>, int> lhs, const pair<pair<int, double>, int> rhs)
	{
		if (lhs.first.first != rhs.first.first)
		{
			return lhs.first.first > rhs.first.first;
		}
		else if (lhs.first.second != rhs.first.second)
		{
			return lhs.first.second < rhs.first.second;
		}
		else
		{
			return (bool)(rand() % 2);
		}
	}
};

struct CustomCompareHJFSQ
{
	bool operator()(const pair<pair<int, double>, int> lhs, const pair<pair<int, double>, int> rhs)
	{
		if ((lhs.first.first / lhs.first.second) != (rhs.first.first / rhs.first.second))
		{
			return (lhs.first.first / lhs.first.second) > (rhs.first.first / rhs.first.second);
		}		
		else if (lhs.first.second != rhs.first.second)
		{
			return lhs.first.second < rhs.first.second;
		}
		else
		{
			return (bool)(rand() % 2);
		}		
	}
};

struct CustomCompareHJFIQ
{
	bool operator()(const pair<pair<int, double>, int> lhs, const pair<pair<int, double>, int> rhs)
	{
		if ((lhs.first.first / lhs.first.second) != (rhs.first.first / rhs.first.second))
		{
			return (lhs.first.first / lhs.first.second) > (rhs.first.first / rhs.first.second);
		}
		else if (lhs.first.second != rhs.first.second)
		{
			return lhs.first.second < rhs.first.second;
		}
		else
		{
			return (bool)(rand() % 2);
		}
	}
};

class sort_indices
{
private:
	double* m_parr;
public:
	sort_indices(double* arr) : m_parr(arr) {};
	bool operator()(int i, int j) const { return m_parr[i] < m_parr[j]; }
};

#endif 

