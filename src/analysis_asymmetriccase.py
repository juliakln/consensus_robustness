"""
Robustness analysis of reaching a stable consensus in the asymmetric case: qx > qy

Functions:
    Analyse model with 1 setting (analyse_model)
    Compare model with different rates qx, qy (compmodels)
    Compare model with different group sizes (compgroupsize)

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

# define file location for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"

# samples for Monte Carlo
samples = 4239

# dictionary of rates with given names to save files accordingly, [qx, qy]
dict_rates = {
    "symm": [1.0, 1.0],
    "smaller": [1.01, 0.99],
    "medium": [1.05, 0.95],
    "medium_2": [1.1, 0.9],
    "larger": [1.4, 0.6]
}

# dictionary of disruptive individuals to save files accordingly
dict_disruptives = {
    "z": "zealots",
    "c": "contrarians"
}

# different linestyles for plotting
linetypes = ['-.', '--', '-', ':', (0, (3, 1, 1, 1)), (0, (3, 10, 1, 10))]

# parameters for font sizes in plots
plt.rcParams.update({
    'legend.fontsize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.titlesize': 18
})


"""
Analyse probability to reach consensus for x and y in 1 setting
Input:
    model: 'voter', 'crossinh' - defines which model of decision-making to use
    N: group size
    stubborn: 'z', 'c' - defines if zealots or contrarians are in the group
    rates: name of rate pair, see dictionary
Output:
    plot of probability curve for X and Y saved in figures/
    elapsed_time: time of the analysis
"""
def analyse_model(model = 'voter', N = 100, stubborn = 'z', rates = 'medium'):
    # measure time of analysis
    start_time = time.time()
    filename = "probs_" + model + "_asym_stableconsensus_" + rates + "rates_" + stubborn + "_" + str(N)
    # compute range of disruptives proportions for which the probabilities are to be calculated 
    stubborn_range = compute_range(N)

    # define distance values based on group size
    distance_base = int(N/10)
    
    # compute probabilities to reach a stable consensus for X and Y
    if model == "voter":
        result_x = stableconsensus_voter_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])
        result_y = stableconsensus_voter_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])
    elif model == "crossinh":
        result_x = stableconsensus_ci_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])
        result_y = stableconsensus_ci_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    fig = plt.figure(figsize=(6,6))
    p1 = [100 * value / N for value in result_x.keys()]
    plt.plot(p1, result_x.values(), 'r', linewidth = 1.5, label = 'X')
    plt.plot(p1, result_y.values(), 'b', linewidth = 1.5, label = 'Y')
    plt.xlim(0,70)
    plt.ylim(0,1)
    plt.xlabel("Percentage of " + dict_disruptives[stubborn] + " in the total group")
    plt.ylabel("Probability to reach consensus")
    plt.title(model + " model, N=" + str(N))
    plt.legend()
    fig.savefig('../figures/' + filename + '.png')
    plt.close()  

    return elapsed_time



"""
Comparison of different rates qx, qy; probabilitiy to reach consensus for X and Y
Input:
    model: 'voter', 'crossinh' - defines which model of decision-making to use
    N: group size
    stubborn: 'z', 'c' - defines if zealots or contrarians are in the group
Output:
    plot of probability curves for X and Y saved in figures/
"""
def compare_rates(model = 'voter', N = 100, stubborn = 'z'):
    # compute range of disruptives proportions for which the probabilities are to be calculated 
    stubborn_range = compute_range(N)

    # define distance values based on group size
    distance_base = int(N/10)

    # save results for each rate pair
    results_x = []
    results_y = []
    labels = [
        rf'$q_X = {rates[0]:.2f},\ q_Y = {rates[1]:.2f}$'
        for rates in dict_rates.values()
    ]
    name = "probs_" + model + "_asym_stableconsensus_"

    if model == "voter":
        # iterate over all rate pairs
        for rates in dict_rates:
            filename = name + rates + "rates_" + stubborn + "_" + str(N)
            results_x.append(stableconsensus_voter_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1]))
            results_y.append(stableconsensus_voter_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1]))

    elif model == "crossinh":
        for rates in dict_rates:
            filename = name + rates + "rates_" + stubborn + "_" + str(N)
            results_x.append(stableconsensus_ci_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1]))
            results_y.append(stableconsensus_ci_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1]))


    fig = plt.figure(figsize=(6.5,6.5))
    for i in range(0, len(results_x)):
        p1 = [100 * value / N for value in results_x[i].keys()]
        plt.plot(p1, results_x[i].values(), linestyle = linetypes[i], color = 'r', label = labels[i])

    for i in range(1, len(results_y)):
        p1 = [100 * value / N for value in results_y[i].keys()]
        plt.plot(p1, results_y[i].values(), linestyle = linetypes[i], color = 'b')

    # dummy_lines = [plt.Line2D([0], [0], color='black', linestyle=style) for style in linetypes]
    # dummy_lines_colors = [plt.Line2D([0], [0], color='red', linestyle='-'), 
    #                   plt.Line2D([0], [0], color='blue', linestyle='-')]
    # legend_labels_colors = ['X', 'Y']  # Labels for color legend


    plt.xlabel('Percentage of ' + dict_disruptives[stubborn] + ' in the group')
    plt.ylabel("Probability to reach consensus")
    plt.xticks(ticks=range(0, 81, 10))

    # legend1 = plt.legend(dummy_lines_colors, legend_labels_colors, 
    #        loc="upper right", bbox_to_anchor=(1, 1))  # Adjust bbox_to_anchor to stack legends

    # plt.gca().add_artist(legend1)  # Add the first legend separately to avoid overwrite

    # plt.legend(dummy_lines, labels, loc="upper right", bbox_to_anchor=(1, 0.82))

    fig.savefig('../figures/stable_' + model + '_asym_' + stubborn + '_' + str(N) + 'font3.png')
    plt.close()   



"""
Comparison of different group sizes N; probabilitiy to reach consensus for X and Y
Input:
    model: 'voter', 'crossinh' - defines which model of decision-making to use
    Ns: array of group size values
    stubborn: 'z', 'c' - defines if zealots or contrarians are in the group
    rates: name of rate pair for qx, qy
Output:
    plot of probability curves for X and Y saved in figures/
"""
def compare_groupsize(model = 'voter', Ns = [100], stubborn = 'z', rates = 'medium'):

    fig = plt.figure(figsize=(6.5,6.5))
    i = 0

    for N in Ns:
        filename = "probs_" + model + "_asym_stableconsensus_" + rates + "rates_" + stubborn + "_" + str(N)
        stubborn_range = compute_range(N)

        # define distance values
        distance_base = int(N/10)

        if model == 'voter':
            result_x = stableconsensus_voter_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])    
            result_y = stableconsensus_voter_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])    
        elif model == 'crossinh':
            result_x = stableconsensus_ci_asym_x(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])    
            result_y = stableconsensus_ci_asym_y(stubborn, N, majority = 50, distance = distance_base, transient = 35, holding = 40, range = stubborn_range, filename = filename, samples = samples, ratex = dict_rates[rates][0], ratey = dict_rates[rates][1])    

        p1 = [100 * value / N for value in result_x.keys()]

        plt.plot(p1, result_x.values(), color = 'r', linestyle = linetypes[i])
        plt.plot(p1, result_y.values(), color = 'b', linestyle = linetypes[i])

        i+=1

    plt.xlabel('Percentage of ' + dict_disruptives[stubborn] + ' in the group')
    plt.ylabel("Probability to reach consensus")

    dummy_lines = [plt.Line2D([0], [0], color='black', linestyle=style) for style in linetypes]
    dummy_lines_colors = [plt.Line2D([0], [0], color='red', linestyle='-'), 
                      plt.Line2D([0], [0], color='blue', linestyle='-')]
    legend_labels_colors = ['X', 'Y']  # Labels for color legend

    legend1 = plt.legend(dummy_lines_colors, legend_labels_colors,
           loc="upper right", bbox_to_anchor=(1, 1))  # Adjust bbox_to_anchor to stack legends

    plt.gca().add_artist(legend1)  # Add the first legend separately to avoid overwrite
    #labels = [str(groupsize) for groupsize in N]
    labels = [rf'$N = {group}$' for group in Ns]
    plt.legend(dummy_lines, labels, loc="upper right", bbox_to_anchor=(1, 0.82))
    
    fig.savefig('../figures/stable_' + model + '_asym_' + stubborn + '_groupsizes_3font.png')
    plt.close()  







def main():
    #analyse_model(model = 'voter', N = 2000, stubborn = 'z', rates = 'medium')

    compare_rates(model = 'voter', N = 100, stubborn = 'z')

    #compare_groupsize(model = 'voter', Ns = [20,50,100,1000,2000,4000], stubborn = 'c', rates = 'medium')



if __name__ == "__main__":
    sys.exit(main())