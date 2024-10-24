"""
Robustness analysis of reaching a stable consensus
- Compare baseline setting for group sizes of 100, 300, 500, 700 and 1000 individuals
- For models with zealots or contrarians

Results are saved in figures/
"""

import os
import sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import time

from utils import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"

colours_groups = ['k', 'g', 'r', 'b', 'c', 'm', 'olive', 'orange', 'darkgreen', 'brown', 'darkviolet']


"""
Compute satisfaction probability of reaching a stable consensus for different consensus settings
    N: total group size
    stubborn: z (zealots), c (contrarians)

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_stable_groups(stubborn = 'z'):

    group_sizes = [10,30,50,70,100,300,500,700,1000]

    start_time = time.time()

#    result = []

    fig = plt.figure(figsize=(6,6))

    i = 0
    for N in group_sizes:

        filename = "probs_stableconsensus2_" + stubborn + "_" + str(N)
        # # 0 to 70% of population are stubborns -> compute 1/2 to have number for each opinion
        # max_stubborn = int(1/2 * (N * 0.7))
        # # test min max_stubborn values and at most 100 values
        # number_values = np.min([max_stubborn+1, 100])
        # # take values between 0 and max value
        # stubborn_range = np.linspace(0, max_stubborn, number_values, dtype=int)
        # #percentages = (((2*stubborn_range) / N) * 100).astype(int)
        
        # We compute the probabilities for maximum 70% of stubborn individuals
        max_stubborn = int(N * 0.7)

        # We want to have around 120 values of zealots for each group size, if possible
        desired_count = 120 
        step = 2
        # Calculate the step size to try to get close to `desired_count` values
        max_possible_values = (max_stubborn - 0) // step + 1
        while max_possible_values > desired_count:
            step += 2
            max_possible_values = (max_stubborn - 0) // step + 1
        # Generate the even numbers within the range
        stubborn_range = np.arange(0, max_stubborn + 1, step)

        # samples for Monte Carlo
        samples = 4239

        # define distance values
        distance_base = int(N/10)

        result = stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples)    
        #perc = np.linspace(0,70,len(result.keys()))
        p1 = [100 * value / N for value in result.keys()]

        plt.plot(p1, result.values(), color = colours_groups[i], label = str(N))
    #    plt.plot(result.keys(), result.values(), color = colours_groups[i], label = str(N))
        i+=1

    #labels = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    #plot_results_groups(stubborn, result, group_sizes, percentages, "Probability to reach a stable consensus", "probs_stableconsensus2_groupsizes")

    plt.xlabel('Amount of contrarians as % of the total group')
    plt.ylabel("Probability to reach a stable consensus")
    plt.legend()
    fig.savefig('../figures/' + "probs_stableconsensus2_c_groupsizes" + '.png')
    plt.close()   

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time



def analyse_switch_groups(stubborn = 'z'):

    group_sizes = [10,30,50,70,100,300,500,700,1000]

    start_time = time.time()

#    result = []

    fig = plt.figure(figsize=(6,6))

    i = 0
    for N in group_sizes:

        filename = "probs_switchconsensus2_" + stubborn + "_" + str(N)
        # # 0 to 70% of population are stubborns -> compute 1/2 to have number for each opinion
        # max_stubborn = int(1/2 * (N * 0.7))
        # # test min max_stubborn values and at most 100 values
        # number_values = np.min([max_stubborn+1, 100])
        # # take values between 0 and max value
        # stubborn_range = np.linspace(0, max_stubborn, number_values, dtype=int)
        # #percentages = (((2*stubborn_range) / N) * 100).astype(int)
        
        # We compute the probabilities for maximum 70% of stubborn individuals
        max_stubborn = int(N)

        # We want to have around 120 values of zealots for each group size, if possible
        desired_count = 120 
        step = 2
        # Calculate the step size to try to get close to `desired_count` values
        max_possible_values = (max_stubborn - 2) // step + 1
        while max_possible_values > desired_count:
            step += 2
            max_possible_values = (max_stubborn - 2) // step + 1
        # Generate the even numbers within the range
        stubborn_range = np.arange(2, max_stubborn + 1, step)

        # samples for Monte Carlo
        samples = 4239

        # define distance values
        distance_base = int(N/10)

        result = switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples)    
        #perc = np.linspace(0,70,len(result.keys()))
        p1 = [100 * value / N for value in result.keys()]

        plt.plot(p1, result.values(), color = colours_groups[i], label = str(N))
    #    plt.plot(result.keys(), result.values(), color = colours_groups[i], label = str(N))
        i+=1

    #labels = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    #plot_results_groups(stubborn, result, group_sizes, percentages, "Probability to reach a stable consensus", "probs_stableconsensus2_groupsizes")

    plt.xlabel('Amount of contrarians as % of the total group')
    plt.ylabel("Probability to switch consensus")
    plt.legend()
    fig.savefig('../figures/' + "probs_switchconsensus2_c_groupsizes" + '.png')
    plt.close()   

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time




"""
Compute satisfaction probability of reaching a stable consensus for different consensus settings
    N: total group size
    stubborn: z (zealots), c (contrarians)

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_stable_groups_voter(stubborn = 'z'):

    group_sizes = [10,30,50,70,100,300,500,700,1000]

    start_time = time.time()

    fig = plt.figure(figsize=(6,6))

    i = 0
    for N in group_sizes:

        filename = "probs_stableconsensus2_voter_" + stubborn + "_" + str(N)
        
        # We compute the probabilities for maximum 70% of stubborn individuals
        max_stubborn = int(N * 0.7)

        # We want to have around 120 values of zealots for each group size, if possible
        desired_count = 120 
        step = 2
        # Calculate the step size to try to get close to `desired_count` values
        max_possible_values = (max_stubborn - 0) // step + 1
        while max_possible_values > desired_count:
            step += 2
            max_possible_values = (max_stubborn - 0) // step + 1
        # Generate the even numbers within the range
        stubborn_range = np.arange(0, max_stubborn + 1, step)

        # samples for Monte Carlo
        samples = 4239

        # define distance values
        distance_base = int(N/10)

        result = stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples)    
        p1 = [100 * value / N for value in result.keys()]

        plt.plot(p1, result.values(), color = colours_groups[i], label = str(N))
        i+=1


    plt.xlabel('Amount of zealots as % of the total group')
    plt.ylabel("Probability to reach a stable consensus")
    plt.legend()
    fig.savefig('../figures/' + "probs_stableconsensus2_voter_z_groupsizes" + '.png')
    plt.close()   

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time




def main():
    #analyse_stable_groups('z')
    #analyse_stable_groups('c')
    #analyse_switch_groups('z')
    #analyse_switch_groups('c')

    analyse_stable_groups_voter('z')


if __name__ == "__main__":
    sys.exit(main())