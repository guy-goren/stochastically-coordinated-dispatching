#!/usr/bin/env python

import os 
import time
import copy
import numpy as np
import multiprocessing
from subprocess import run
from pathlib import Path

###############################################################################
###############################################################################

def worker(simArgs, simType, sync=None, d=None):
                           
    command = simArgs['path'] + ' '\
            + str(simArgs['number_of_servers']) + ' '\
            + str(simArgs['number_of_dispatchers']) + ' '\
            + str(simArgs['sim_time']) + ' '\
    
    for l in simArgs['arrival_parameters']:
        
        command += ' '
        command += str(l)
        
    for s in simArgs['departure_parameters']:
        
        command += ' '
        command += str(s)
    
    command += ' ' 
    command += simType
    
    if simType == 'lsq':

        command += ' ' 
        command += str(d)
         
    print(command)
    
    start = time.time()
    run(command)	
    end = time.time()
    
    runType = end - start
               
    return runType

###############################################################################
###############################################################################

loads = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99]
   
###############################################################################
###############################################################################

def complete_state_information(args):
    
    func_args = copy.deepcopy(args)
    
    totalDepartureRate = sum([s for s in args['departure_parameters']]) 
    
    for load in loads:
        
        totalArrivalRate = totalDepartureRate*load
           
        func_args['arrival_parameters'] = [float(totalArrivalRate)/func_args['number_of_dispatchers']]*func_args['number_of_dispatchers']
        
        print(worker(func_args, 'complete_information'))
            
###############################################################################
###############################################################################

def weighted_random(args):
    
    func_args = copy.deepcopy(args)
    
    totalDepartureRate = sum([s for s in args['departure_parameters']]) 
    
    for load in loads:
        
        totalArrivalRate = totalDepartureRate*load
           
        func_args['arrival_parameters'] = [float(totalArrivalRate)/func_args['number_of_dispatchers']]*func_args['number_of_dispatchers']
        
        print(worker(func_args, 'wr'))
            
###############################################################################
############################################################################### 

def jiq(args):
    
    func_args = copy.deepcopy(args)
    
    totalDepartureRate = sum([s for s in args['departure_parameters']]) 
    
    for load in loads:
        
        totalArrivalRate = totalDepartureRate*load
           
        func_args['arrival_parameters'] = [float(totalArrivalRate)/func_args['number_of_dispatchers']]*func_args['number_of_dispatchers']
                                
        print(worker(func_args, 'jiq'))
         
###############################################################################
############################################################################### 

def lsq(args):
    
    func_args = copy.deepcopy(args)
    
    totalDepartureRate = sum([s for s in args['departure_parameters']]) 
    
    for load in loads:
        
        totalArrivalRate = totalDepartureRate*load
           
        func_args['arrival_parameters'] = [float(totalArrivalRate)/func_args['number_of_dispatchers']]*func_args['number_of_dispatchers']
                                
        print(worker(func_args, 'lsq', d=2))
            
###############################################################################
############################################################################### 
            
if __name__ == "__main__":

    ###########################################################################
    ###########################################################################
    
    read_service_rates_from_file = False
    
    heterogeneity_level = 100
    
    ###########################################################################
    ###########################################################################
    
    simTypes = ['complete_state_information', 'jiq', 'lsq']
    
    ###########################################################################
    ###########################################################################
    
    for number_of_servers, number_of_dispatchers in [(300,10), (400,10)]:#[(100,5), (100, 10), (200, 10), (200,20)]:
        
    
        args = {}
        
        
        dir_path = os.path.dirname(os.path.realpath(__file__))
        
        args['path'] = str(Path(dir_path).parent) + "\\x64\Release\HeterogeneousLoadBalancing.exe"
        
        args['number_of_servers'] = number_of_servers    
        args['number_of_dispatchers'] = number_of_dispatchers
    
        args['sim_time'] = 1000#00
        
        '''
        departure_parameters_1 = [2.0]*int(args['number_of_servers']*0.1) 
        departure_parameters_2 = [20.0]*int(args['number_of_servers']*0.9) 
        
        args['departure_parameters'] = departure_parameters_1 + departure_parameters_2
        '''
        
        if read_service_rates_from_file:
            
            args['departure_parameters'] = []
            
            path = Path(os.getcwd())
            
            filepath = str(path.parent) + "\\WindowsEnv\\random_1_{}\\".format(heterogeneity_level) + "service rates_n_{}_m_{}.txt".format(number_of_servers, number_of_dispatchers)
            with open(filepath) as fp:
                lines = fp.readlines()
                for line in lines:
                    args['departure_parameters'].append(float(line.split()[0]))
            print("loaded service rates: ", args['departure_parameters'])

        else:
            
            args['departure_parameters'] = [1 + (heterogeneity_level - 1) * np.random.random() for _ in range(args['number_of_servers'])]
        
        ######################################################################
        
        f = open("service rates_n_{}_m_{}.txt".format(number_of_servers, number_of_dispatchers), "a")
        for i in range(args['number_of_servers']):
            f.write("{}\n".format(args['departure_parameters'][i]))
        f.close()

        ######################################################################
                    
        procs = []
                                    
        procs.append(multiprocessing.Process(target=complete_state_information, args=(args,)))
        procs.append(multiprocessing.Process(target=weighted_random, args=(args,)))
        procs.append(multiprocessing.Process(target=jiq, args=(args,)))
        procs.append(multiprocessing.Process(target=lsq, args=(args,)))
                    
        for p in procs:
            p.start()
        
        
        for p in procs:
            p.join()    
        
