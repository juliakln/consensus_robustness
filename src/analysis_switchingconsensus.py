"""
Robustness analysis of switching consensus
- For group sizes of 100, 500, and 1000 individuals
- For models with zealots, contrarians, or both

Results are saved in figures/
"""

import os
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


"""
Compute satisfaction probability of switching consensus for different consensus settings
    N: total group size
    stubborn: z (zealots), c (contrarians)

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_switch(N = 100, stubborn = 'z'):
    start_time = time.time()

    filename = "probs_switchconsensus2_" + stubborn + "_" + str(N)
 
     # We compute the probabilities for maximum 100% of stubborn individuals
    max_stubborn = int(N)
    #max_stubborn = int(1/2 * N)

    # We want to have around 120 values of zealots for each group size, if possible
    desired_count = 120 

    step = 2
    # Calculate the step size to try to get close to `desired_count` values
    max_possible_values = (max_stubborn - 2) // step + 1
    while max_possible_values > desired_count:
        step += 2
        max_possible_values = (max_stubborn - 2) // step + 1

    # if max_possible_values <= desired_count:
    #     step = 2
    # else:
    #     step = 2 * ((max_stubborn - 0) // (2 * desired_count))
    
    # Generate the even numbers within the range
    stubborn_range = np.arange(2, max_stubborn + 1, step)

    # # test min max_stubborn values and at most 100 values
    # number_values = np.min([max_stubborn+1, 100])
    # # take values between 1 and max value
    # stubborn_range = np.linspace(1, max_stubborn, number_values, dtype=int)
    # percentages = (((2*stubborn_range) / N) * 100).astype(int)
    # samples for Monte Carlo
    samples = 4239

    # define distance values
    distance_base = int(N/10) # 5, 7, 12, 15; base: 10
    low1  = int(N*0.05)
    low2  = int(N*0.07)
    high1 = int(N*0.12)
    high2 = int(N*0.15)

    result = []
    result.append(switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = low1, transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = low2, transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = high1,  transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = high2, transient = 35, holding = 10, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = low1, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = low2, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = high1, range = stubborn_range, filename = filename, samples = samples))
    result.append(switchconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = high2, range = stubborn_range, filename = filename, samples = samples))

    labels = ['Baseline', 'd='+str(low1), 'd='+str(low2), 'd='+str(high1), 'd='+str(high2),
                          's='+str(low1), 's='+str(low2), 's='+str(high1), 's='+str(high2)]

    plot_results_switch(stubborn, result, labels, N, "Probability to switch consensus", filename)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time



# def analyse_switch_100(stubborn):

#     N = 100
#     filename = "probs_switchconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1, 50, 50)
#     samples = 4239

#     result = []
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 35, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 65, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 1,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 20, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=1', 'd=20', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability of switching consensus", filename)


# def analyse_switch_250(stubborn):

#     N = 250
#     filename = "probs_switchconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1, 98, 50)
#     samples = 4239

#     result = []
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 35, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 65, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 2,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 25, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 25, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=2', 'd=50', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability of switching consensus", filename)


# def analyse_switch_500(stubborn):

#     N = 500
#     filename = "probs_switchconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1,244,49)
#     samples = 1060

#     result = []
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 35, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 65, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 5,   transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50,  transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50,  transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=5', 'd=100', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability of switching consensus", filename)


# def analyse_switch_1000(stubborn):

#     N = 1000
#     filename = "probs_switchconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1,491,50)
#     samples = 1060

#     result = []
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 35, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 65, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 10,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 100, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 100, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(switchconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=10', 'd=200', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability of switching consensus", filename)



"""
Compute satisfaction probability of switching consensus for baseline consensus settings

    Output: contour plot over #stubborn individuals 
"""
def analyse_switch_100_both():

    N = 100
    filename = "probs_switchconsensus_both_" + str(N)
    zealots = np.linspace(1,49,25)
    contrarians = np.linspace(1,19,10)
    samples = 1060

    result = switchconsensus_both(N, majority = 50, distance = 10, transient = 35, holding = 10, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)


def analyse_switch_500_both():

    N = 500
    filename = "probs_switchconsensus_both_" + str(N)
    zealots = np.linspace(1,241,25)
    contrarians = np.linspace(1,100,10)
    samples = 1060

    result = switchconsensus_both(N, majority = 50, distance = 50, transient = 35, holding = 10, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)


def analyse_switch_1000_both():

    N = 1000
    filename = "probs_switchconsensus_both_" + str(N)
    zealots = np.linspace(1,481,25)
    contrarians = np.linspace(1,496,10)
    samples = 1060

    result = switchconsensus_both(N, majority = 50, distance = 100, transient = 35, holding = 10, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)