#ifndef SIMULATIONS_
#define SIMULATIONS_

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void complete_information(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType);

void complete_information_TWF(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void classical_jiq(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, string dispatcherType);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void lsq_sample(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, int d, string dispatcherType);

void lsq_sample_het(int n, int m, int time, double* arrival_parameters, double* departure_parameters, int seed, int fn_resolution, int d, string dispatcherType);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

#endif