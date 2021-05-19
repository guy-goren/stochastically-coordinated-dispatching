#include "server.h"

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

Server::Server(int id) :
	m_serverID(id)
{
}

Server::~Server()
{
}

int Server::getQueueSize()
{
	return (int)m_queue.size();
}

void Server::insertJobs(int jobs, int currentTime)
{
	for (int j = 0; j < jobs; ++j)
	{
		job_t job;
		job.creation_time = currentTime;
		m_queue.push(job);
	}
}

vector<int> Server::serveJobsAndReturnCreationTimes(int jobs)
{
	int toComplete = (int)m_queue.size() < jobs ? (int)m_queue.size() : jobs;

	vector<int> creationTimes;
	for (int j = 0; j < toComplete; ++j)
	{
		int creationTime = m_queue.front().creation_time;
		creationTimes.push_back(creationTime);
		m_queue.pop();
	}	
	return creationTimes;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

ServerBasic::ServerBasic(int id, int n) :
	Server(id),
	m_NumberOfServers(n)
{
	m_localState = new SizeTimestampPair[n]();
	for (int i = 0; i < n; ++i)
	{
		m_localState[i].first = 0;
		m_localState[i].second = NULL_TIMESTAMP;
	}
}

ServerBasic::~ServerBasic()
{
	delete[] m_localState;
}

void ServerBasic::getInformation(SizeTimestampPair * sizeTimestampPairs, int currentTime)
{
	m_localState[m_serverID].first = (int)m_queue.size();
	m_localState[m_serverID].second = currentTime;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		sizeTimestampPairs[i] = m_localState[i];
	}
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

ServerAdvanced::ServerAdvanced(int id, int n) :
	ServerBasic(id, n)
{
}

ServerAdvanced::~ServerAdvanced()
{
}

void ServerAdvanced::setInformation(const SizeTimestampPair * sizeTimestampPairs, int currentTime)
{
	m_localState[m_serverID].first = (int)m_queue.size();
	m_localState[m_serverID].second = currentTime;

	for (int i = 0; i < m_NumberOfServers; ++i)
	{
		// if the local timestamp is older update both queue length and its timestamp
		if (m_localState[i].second < sizeTimestampPairs[i].second)
		{
			m_localState[i] = sizeTimestampPairs[i];
		}
	}
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

ServerJIQ::ServerJIQ(int id) :
	Server(id)
{
	m_idled = true;
}

ServerJIQ::~ServerJIQ()
{
}

vector<int> ServerJIQ::serveJobsAndReturnCreationTimesJIQ(int jobs)
{
	bool goingToWork = jobs + (int)m_queue.size() > 0 ? true : false;
	vector<int> creationTimes = serveJobsAndReturnCreationTimes(jobs);
	if (goingToWork && ((int)m_queue.size() == 0))
	{
		m_idled = true;
	}
	return creationTimes;
}

bool ServerJIQ::gotIdle()
{
	if (m_idled)
	{
		m_idled = false;
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

