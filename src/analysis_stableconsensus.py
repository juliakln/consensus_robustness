"""
Robustness analysis of reaching a stable consensus
- For group sizes of 100, 500, and 1000 individuals
- For models with zealots, contrarians, or both

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


"""
Compute satisfaction probability of reaching a stable consensus for different consensus settings
    N: total group size
    stubborn: z (zealots), c (contrarians)

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_stable(N = 100, stubborn = 'z'):

    start_time = time.time()

    filename = "probs_stableconsensus2_" + stubborn + "_" + str(N)

    # We compute the probabilities for maximum 70% of stubborn individuals
    max_stubborn = int(N * 0.7)

    # We want to have around 120 values of zealots for each group size, if possible
    desired_count = 120 

    # Calculate the maximum possible step size to get at most `desired_count` values
    #max_step = (max_stubborn - 0) // (desired_count - 1)
    
    # # Adjust step size to be even
    # step = 2 * (max_step // 2)
    
    # # If step size is less than 2, set it to 2 to ensure even numbers
    # if step < 2:
    #     step = 2
    
    # # Generate even numbers within [0,(70%N)] -> even, because we want to have symmetric amounts of zealots for both options, Z=Zx+Zy, Zx=Zy
    # stubborn_range = np.arange(0, max_stubborn + 1, step)
    
    # # If we have fewer values than desired, reduce the step size
    # while len(stubborn_range) < desired_count:
    #     step -= 2
    #     stubborn_range = np.arange(0, max_stubborn + 1, step)

    step = 2
    # Calculate the step size to try to get close to `desired_count` values
    max_possible_values = (max_stubborn - 0) // step + 1
    while max_possible_values > desired_count:
        step += 2
        max_possible_values = (max_stubborn - 0) // step + 1

    # if max_possible_values <= desired_count:
    #     step = 2
    # else:
    #     step = 2 * ((max_stubborn - 0) // (2 * desired_count))
    
    # Generate the even numbers within the range
    stubborn_range = np.arange(0, max_stubborn + 1, step)
    
    # If we have fewer values than desired, incrementally reduce the step size
    # while len(stubborn_range) < desired_count and step > 2:
    #     step -= 2
    #     stubborn_range = np.arange(0, max_stubborn + 1, step)


    # # 0 to 70% of population are stubborns -> compute 1/2 to have number for each opinion
    # max_stubborn = int(1/2 * (N * 0.7))
    # # test min max_stubborn values and at most 100 values
    # number_values = np.min([max_stubborn+1, 100])
    # # take values between 0 and max value
    # stubborn_range = np.linspace(0, max_stubborn, number_values, dtype=int)
    # percentages = (((2*stubborn_range) / N) * 100).astype(int)
    
    # samples for Monte Carlo
    samples = 4239

    # define distance values
    distance_base = int(N/10)
    distance_low  = int(N/100)
    distance_high = int(N/5)

    result = []
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 35, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 65, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_low,  transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_high, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 20, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 50, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 25, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 55, range = stubborn_range, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    plot_results(stubborn, result, labels, N, "Probability to reach a stable consensus", filename)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time


def analyse_stable_voter(N = 100, stubborn = 'z'):

    start_time = time.time()

    filename = "probs_voter_stableconsensus2_" + stubborn + "_" + str(N)

    # We compute the probabilities for maximum 70% of stubborn individuals
    max_stubborn = int(N * 0.7)

    # We want to have around 120 values of zealots for each group size, if possible
    desired_count = 120 

    # Calculate the maximum possible step size to get at most `desired_count` values
    #max_step = (max_stubborn - 0) // (desired_count - 1)
    
    # # Adjust step size to be even
    # step = 2 * (max_step // 2)
    
    # # If step size is less than 2, set it to 2 to ensure even numbers
    # if step < 2:
    #     step = 2
    
    # # Generate even numbers within [0,(70%N)] -> even, because we want to have symmetric amounts of zealots for both options, Z=Zx+Zy, Zx=Zy
    # stubborn_range = np.arange(0, max_stubborn + 1, step)
    
    # # If we have fewer values than desired, reduce the step size
    # while len(stubborn_range) < desired_count:
    #     step -= 2
    #     stubborn_range = np.arange(0, max_stubborn + 1, step)

    step = 2
    # Calculate the step size to try to get close to `desired_count` values
    max_possible_values = (max_stubborn - 0) // step + 1
    while max_possible_values > desired_count:
        step += 2
        max_possible_values = (max_stubborn - 0) // step + 1

    # if max_possible_values <= desired_count:
    #     step = 2
    # else:
    #     step = 2 * ((max_stubborn - 0) // (2 * desired_count))
    
    # Generate the even numbers within the range
    stubborn_range = np.arange(0, max_stubborn + 1, step)
    
    # If we have fewer values than desired, incrementally reduce the step size
    # while len(stubborn_range) < desired_count and step > 2:
    #     step -= 2
    #     stubborn_range = np.arange(0, max_stubborn + 1, step)


    # # 0 to 70% of population are stubborns -> compute 1/2 to have number for each opinion
    # max_stubborn = int(1/2 * (N * 0.7))
    # # test min max_stubborn values and at most 100 values
    # number_values = np.min([max_stubborn+1, 100])
    # # take values between 0 and max value
    # stubborn_range = np.linspace(0, max_stubborn, number_values, dtype=int)
    # percentages = (((2*stubborn_range) / N) * 100).astype(int)
    
    # samples for Monte Carlo
    samples = 4239

    # define distance values
    distance_base = int(N/10)
    distance_low  = int(N/100)
    distance_high = int(N/5)

    result = []
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 35, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 65, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_low,  transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_high, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 20, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 50, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 25, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 55, range = stubborn_range, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    plot_results(stubborn, result, labels, N, "Probability to reach a stable consensus", filename)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return elapsed_time

# def analyse_stable_100(stubborn):

#     N = 100
#     filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1, 35, 35)
#     samples = 4239

#     result = []
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 35, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 65, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 1,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 20, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=1', 'd=20', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)


# def analyse_stable_250(stubborn):

#     N = 250
#     filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1, 89, 45)
#     samples = 4239

#     result = []
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 35, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 65, distance = 25, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 2,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 25, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 25, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 25, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=2', 'd=50', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)



# def analyse_stable_500(stubborn):

#     N = 500
#     filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1,148,50)
#     samples = 1060

#     result = []
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 35, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 65, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 5,   transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=5', 'd=100', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)


# def analyse_stable_1000(stubborn):

#     N = 1000
#     filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
#     zealots = np.linspace(1,344,50)
#     samples = 1060

#     result = []
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 35, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 65, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 10,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
#     result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

#     labels = ['Baseline', 'm=35', 'm=65', 'd=10', 'd=200', 't=20', 't=50', 'h=25', 'h=55']

#     plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)



"""
Compute satisfaction probability of reaching a stable consensus for baseline consensus settings

    Output: contour plot over #stubborn individuals 
"""
def analyse_stable_100_both():

    N = 100
    filename = "probs_stableconsensus_both_" + str(N)
    zealots = np.linspace(1,49,25)
    contrarians = np.linspace(1,19,10)
    samples = 1060

    result = stableconsensus_both(N, majority = 50, distance = 10, transient = 35, holding = 40, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)


def analyse_stable_500_both():

    N = 500
    filename = "probs_stableconsensus_both_" + str(N)
    zealots = np.linspace(1,241,25)
    contrarians = np.linspace(1,100,10)
    samples = 1060

    result = stableconsensus_both(N, majority = 50, distance = 50, transient = 35, holding = 40, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)


def analyse_stable_1000_both():

    N = 1000
    filename = "probs_stableconsensus_both_" + str(N)
    zealots = np.linspace(1,481,25)
    contrarians = np.linspace(1,496,10)
    samples = 1060

    result = stableconsensus_both(N, majority = 50, distance = 100, transient = 35, holding = 40, range_z = zealots, range_c = contrarians, filename = filename, samples = samples)

    plot_results_2dim(result, zealots, contrarians, filename)




def analyse_stable_CIandVoter(N = 100):

    filename = "probs_stableconsensus2_z_" + str(N)
    stubborn = 'z'
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
    distance_low  = int(N/100)
    distance_high = int(N/5)

    result = []
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 35, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 65, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_low,  transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_high, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 20, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 50, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 25, range = stubborn_range, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 55, range = stubborn_range, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    #plot_results(stubborn, result, labels, N, "Probability to reach a stable consensus", filename)


    filename = "probs_voter_stableconsensus2_" + stubborn + "_" + str(N)

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
    distance_low  = int(N/100)
    distance_high = int(N/5)

    result2 = []
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 35, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 65, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_low,  transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_high, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 20, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 50, holding = 40, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 25, range = stubborn_range, filename = filename, samples = samples))
    result2.append(stableconsensus_voter(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 55, range = stubborn_range, filename = filename, samples = samples))

    labels2 = ['Baseline', 'm=35', 'm=65', 'd='+str(distance_low), 'd='+str(distance_high), 't=20', 't=50', 'h=25', 'h=55']

    #plot_results(stubborn, result, labels, N, "Probability to reach a stable consensus", filename)

    fig = plt.figure(figsize=(6,6))
    for i in range(1, len(result)):
#        perc = np.linspace(0,70,len(results[i].keys()))
        p1 = [100 * value / N for value in result[i].keys()]
        plt.plot(p1, result[i].values(), linestyle = linetypes[i], color = colours[i], label = labels[i])
    #perc = np.linspace(0,70,len(results[0].keys()))
    p1 = [100 * value / N for value in result[0].keys()]
    plt.plot(p1, result[0].values(), 'k', linewidth = 1.5, label = 'Baseline')

    for i in range(1, len(result2)):
#        perc = np.linspace(0,70,len(results[i].keys()))
        p1 = [100 * value / N for value in result2[i].keys()]
        plt.plot(p1, result2[i].values(), linestyle = linetypes[i], color = colours[i], alpha=0.5)
    #perc = np.linspace(0,70,len(results[0].keys()))
    p1 = [100 * value / N for value in result[0].keys()]
    plt.plot(p1, result2[0].values(), 'k', linewidth = 1.5, alpha=0.5)

    plt.xlabel('Amount of zealots as % of the total group')
    plt.ylabel("Probability to reach a stable consensus")
    plt.legend()
    fig.savefig('../figures/stable_voter_crossinhibition_z_100.png')
    plt.close()   

    return 3







def main():
    #analyse_stable_voter(N=100, stubborn='z')
    analyse_stable_CIandVoter()

if __name__ == "__main__":
    sys.exit(main())