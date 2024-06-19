"""
Robustness analysis of reaching a stable consensus
- For group sizes of 100, 500, and 1000 individuals
- For models with zealots, contrarians, or both

Results are saved in figures/
"""

import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

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
    stubborn: z (zealots), c (contrarians)

    Output: lineplot over #stubborn individuals for different settings
"""
def analyse_stable_100(stubborn):

    N = 100
    filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
    zealots = np.linspace(1, 35, 35)
    samples = 4239

    result = []
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 35, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 65, distance = 10, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 1,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 20, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd=1', 'd=20', 't=20', 't=50', 'h=25', 'h=55']

    plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)



def analyse_stable_500(stubborn):

    N = 500
    filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
    zealots = np.linspace(1,148,50)
    samples = 1060

    result = []
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 35, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 65, distance = 50,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 5,   transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 50,  transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd=5', 'd=100', 't=20', 't=50', 'h=25', 'h=55']

    plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)


def analyse_stable_1000(stubborn):

    N = 1000
    filename = "probs_stableconsensus_" + stubborn + "_" + str(N)
    zealots = np.linspace(1,344,50)
    samples = 1060

    result = []
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 35, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 65, distance = 100, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 10,  transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 20, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 100, transient = 50, holding = 40, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 25, range = zealots, filename = filename, samples = samples))
    result.append(stableconsensus(stubborn, N, majority = 50, distance = 200, transient = 35, holding = 55, range = zealots, filename = filename, samples = samples))

    labels = ['Baseline', 'm=35', 'm=65', 'd=10', 'd=200', 't=20', 't=50', 'h=25', 'h=55']

    plot_results(stubborn, result, labels, zealots, "Probability to reach a stable consensus", filename)



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
