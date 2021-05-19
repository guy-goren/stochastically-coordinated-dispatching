#include "defs.h"
#include "logger.h"

Logger::Logger(int simulation_time, int sample_resolution) : st_(simulation_time), sr_(sample_resolution)
{
	job_jcts_cumsum_ = 0;
	completed_jobs_ = 0;
	queue_lengths_cumsum_ = 0;

	nsu_ = simulation_time / sample_resolution;
	ta_queue_lengths_.push_back(0);
}

Logger::~Logger()
{
}

void Logger::updateMQL(int * queueLengths, int n, int time)
{
	for (int i = 0; i < n; ++i)
	{
		queue_lengths_cumsum_ += queueLengths[i];
	}

	if (nsu_ > 0)
	{
		--nsu_;
	}
	else
	{
		ta_queue_lengths_.push_back(queue_lengths_cumsum_ / (long double)time);
		nsu_ = st_ / sr_;
	}
}

void Logger::updateJCT(int job_creation_time, int time)
{
	++completed_jobs_;

	int jct = time - job_creation_time + 1;

	job_jcts_cumsum_ += jct;
	jct_map_[jct] += 1;
}

void Logger::updateGDT(long long getDestinations_nanosec_exectime)
{
	gdt_map_[getDestinations_nanosec_exectime] += 1;
}

void Logger::printStatus(int time, int ps_resolution, int current_sum)
{
	// Collect stats and print
	if (time % ps_resolution == 0) {
		cout << "Time: " << time << " Jobs: " << queue_lengths_cumsum_ / (long double)time << " Total jobs " << completed_jobs_ + current_sum << endl;
	}
}

void Logger::writeResultsToFiles(string fn, int time)
{
	ofstream data_file;
	string data_file_fn;
	data_file_fn.append(fn);
	data_file_fn.append("_mean_queue_lengths.txt");
	data_file.open(data_file_fn, std::ofstream::out | std::ofstream::app);
	data_file << queue_lengths_cumsum_ / (long double)time << endl;
	data_file.close();

	ofstream vector_data_file;
	string vector_data_file_fn;
	vector_data_file_fn.append(fn);
	vector_data_file_fn.append("_mean_queue_lengths_evolution.txt");
	vector_data_file.open(vector_data_file_fn, std::ofstream::out | std::ofstream::app);
	for (vector<long double>::iterator it = ta_queue_lengths_.begin(); it != ta_queue_lengths_.end(); ++it) {
		vector_data_file << *it << endl;
	}
	vector_data_file.close();

	ofstream jct_data_file;
	string jct_data_file_fn;
	jct_data_file_fn.append(fn);
	jct_data_file_fn.append("_mean_jcts.txt");
	jct_data_file.open(jct_data_file_fn, ofstream::out | ofstream::app);
	if (completed_jobs_ > 0)
	{
		jct_data_file << job_jcts_cumsum_ / completed_jobs_ << endl;
	}
	jct_data_file.close();

	ofstream jct_map_data_file;
	string jct_map_data_file_fn;
	jct_map_data_file_fn.append(fn);
	jct_map_data_file_fn.append("_jct_histogram.txt");
	jct_map_data_file.open(jct_map_data_file_fn, ofstream::out | ofstream::app);
	for (auto& t : jct_map_)
	{
		jct_map_data_file << t.first << "\t" << t.second << endl;
	}
	jct_map_data_file.close();

	ofstream gdt_map_data_file;
	string gdt_map_data_file_fn;
	gdt_map_data_file_fn.append(fn);
	gdt_map_data_file_fn.append("_gdt_histogram.txt");
	gdt_map_data_file.open(gdt_map_data_file_fn, ofstream::out | ofstream::app);
	for (auto& t : gdt_map_)
	{
		gdt_map_data_file << t.first << "\t" << t.second << endl;
	}
	gdt_map_data_file.close();
}

