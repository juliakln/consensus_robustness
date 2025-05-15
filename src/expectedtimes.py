"""
Expected times analysis of reaching and holding a stable consensus

Functions:
    Verify reaching/holding property with Prism (verify_reaching_prism)
    Read expected times from Prism output files (read_reaching_prism)
    Plot expected times (plot_reaching)

Results are saved in data/times/. and figures/
"""



import os
import sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import time
import subprocess

from utils import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# dictionary of rates with given names to save files accordingly, [qx, qy], for N=100
dict_rates = {
    "symm": [0.01, 0.01],
    "smaller": [0.0101, 0.0099],
    "medium": [0.0105, 0.0095],
    "medium_2": [0.011, 0.009],
    "larger": [0.014, 0.006]
}



def verify_reaching_prism_zealots(rates):

    # amount of zealots from 2 to 90 (z = zx = zy)
    zealots_range = np.arange(1, 46, 1) 

    # compute expected time for reaching consensus with Prism
    for zealots in zealots_range:
        result = "/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_reach/reaching_" + str(zealots) + ".txt"
        if not os.path.exists(result):
            os.makedirs(os.path.dirname(result), exist_ok=True)
            resultfile = open(result ,"w")
            resultfile.close()
            resultfile = open(result, "r+")
            prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/consensus_robustness/models/tanja_cross_inhibition.pm /Users/juliaklein/Documents/consensus_robustness/models/tanja_CI.pctl -const Zx=" + str(zealots) + ",Zy=" + str(zealots) + ",qx=" + str(dict_rates[rates][0]) + ",qy=" + str(dict_rates[rates][1]) + " -prop 1 -sparse -maxiters 1000000 -exportresults " + result
            prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
            resultfile.close()


def verify_reaching_prism_contrarians(rates):

    # amount of contrarians from 2 to 90 (z = zx = zy)
    contrarians_range = np.arange(1, 46, 1) 

    # compute expected time for reaching consensus with Prism
    for contrarians in contrarians_range:
        result = "/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_reach_contrarians/reaching_" + str(contrarians) + ".txt"
        if not os.path.exists(result):
            os.makedirs(os.path.dirname(result), exist_ok=True)
            resultfile = open(result ,"w")
            resultfile.close()
            resultfile = open(result, "r+")
            prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/consensus_robustness/models/tanja_cross_inhibition_contrarians.pm /Users/juliaklein/Documents/consensus_robustness/models/tanja_CI.pctl -const cHalf=" + str(contrarians) + ",qx=" + str(dict_rates[rates][0]) + ",qy=" + str(dict_rates[rates][1]) + " -prop 1 -sparse -maxiters 1000000 -exportresults " + result
            prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
            resultfile.close()



def verify_holding_prism_zealots(rates):

    # amount of zealots from 2 to 90 (z = zx = zy)
    zealots_range = np.arange(1, 46, 1) 

    # compute expected time for reaching consensus with Prism
    for zealots in zealots_range:
        result = "/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_hold/holding_" + str(zealots) + ".txt"
        if not os.path.exists(result):
            os.makedirs(os.path.dirname(result), exist_ok=True)
            resultfile = open(result ,"w")
            resultfile.close()
            resultfile = open(result, "r+")
            prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/consensus_robustness/models/tanja_cross_inhibition.pm /Users/juliaklein/Documents/consensus_robustness/models/tanja_CI.pctl -const Zx=" + str(zealots) + ",Zy=" + str(zealots) + ",qx=" + str(dict_rates[rates][0]) + ",qy=" + str(dict_rates[rates][1]) + " -prop 2 -sparse -maxiters 10000000 -exportresults " + result
            prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
            resultfile.close()


def verify_holding_prism_contrarians(rates):

    # amount of contrarians from 2 to 90 (z = zx = zy)
    contrarians_range = np.arange(1, 46, 1) 

    # compute expected time for reaching consensus with Prism
    for contrarians in contrarians_range:
        result = "/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_hold_contrarians/holding_" + str(contrarians) + ".txt"
        if not os.path.exists(result):
            os.makedirs(os.path.dirname(result), exist_ok=True)
            resultfile = open(result ,"w")
            resultfile.close()
            resultfile = open(result, "r+")
            prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/consensus_robustness/models/tanja_cross_inhibition_contrarians.pm /Users/juliaklein/Documents/consensus_robustness/models/tanja_CI.pctl -const cHalf=" + str(contrarians) + ",qx=" + str(dict_rates[rates][0]) + ",qy=" + str(dict_rates[rates][1]) + " -prop 2 -sparse -maxiters 100000000 -exportresults " + result
            prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
            resultfile.close()



def read_reaching_prism_zealots(rates):
    # save expected times for each value of p and experiment
    times = defaultdict(list)
    # read outcomes for all parameter values and compute number of satisfactions
    for dirpath, dirs, files in os.walk("/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_reach/"):
        for file in files:
            if file.startswith("reaching"):
                zealots = round(int((file.split("_")[1]).rsplit(".", 1)[0]), 0)
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.readline()
                    for last_line in f:
                        pass
                    try:
                        t = float(last_line)
                    except ValueError:
                        t = np.nan
                    times[zealots] = t

    # Training data
    paramValueSet = []
    paramValueOutputs = []
    for key in sorted(times):
        paramValueSet.append(key)
        paramValueOutputs.append(times[key])

    return paramValueSet, paramValueOutputs



def read_reaching_prism_contrarians(rates):
    # save expected times for each value of p and experiment
    times = defaultdict(list)
    # read outcomes for all parameter values and compute number of satisfactions
    for dirpath, dirs, files in os.walk("/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_reach_contrarians/"):
        for file in files:
            if file.startswith("reaching"):
                zealots = round(int((file.split("_")[1]).rsplit(".", 1)[0]), 0)
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.readline()
                    for last_line in f:
                        pass
                    try:
                        t = float(last_line)
                    except ValueError:
                        t = np.nan
                    times[zealots] = t

    # Training data
    paramValueSet = []
    paramValueOutputs = []
    for key in sorted(times):
        paramValueSet.append(key)
        paramValueOutputs.append(times[key])

    return paramValueSet, paramValueOutputs



def read_holding_prism_zealots(rates):
    # save expected times for each value of p and experiment
    times = defaultdict(list)
    # read outcomes for all parameter values and compute number of satisfactions
    for dirpath, dirs, files in os.walk("/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_hold/"):
        for file in files:
            if file.startswith("holding"):
                zealots = round(int((file.split("_")[1]).rsplit(".", 1)[0]), 0)
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.readline()
                    for last_line in f:
                        pass
                    try:
                        t = float(last_line)
                    except ValueError:
                        t = np.nan
                    times[zealots] = t

    # Training data
    paramValueSet = []
    paramValueOutputs = []
    for key in sorted(times):
        paramValueSet.append(key)
        paramValueOutputs.append(times[key])

    return paramValueSet, paramValueOutputs



def read_holding_prism_contrarians(rates):
    # save expected times for each value of p and experiment
    times = defaultdict(list)
    # read outcomes for all parameter values and compute number of satisfactions
    for dirpath, dirs, files in os.walk("/Users/juliaklein/Documents/consensus_robustness/data/times/" + rates + "_hold_contrarians/"):
        for file in files:
            if file.startswith("holding"):
                zealots = round(int((file.split("_")[1]).rsplit(".", 1)[0]), 0)
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.readline()
                    for last_line in f:
                        pass
                    try:
                        t = float(last_line)
                    except ValueError:
                        t = np.nan
                    times[zealots] = t

    # Training data
    paramValueSet = []
    paramValueOutputs = []
    for key in sorted(times):
        paramValueSet.append(key)
        paramValueOutputs.append(times[key])

    return paramValueSet, paramValueOutputs



def plot_reaching(rates, zealots, reachingtimes_z, contrarians, reachingtimes_c):

    fig = plt.figure(figsize=(6,6))
    p1 = [2 * value for value in zealots]
    plt.plot(p1, np.log(reachingtimes_z), 'b', linewidth = 1.5, label = 'Zealots')
    plt.plot(p1, np.log(reachingtimes_c), 'r', linewidth = 1.5, label = 'Contrarians')
    plt.xlim(0,100)
    plt.ylim(-3, 8)
    plt.xlabel("Amount of disruptive individuals")
    plt.ylabel("Expected time to reach consensus in log scale")
    plt.title("Asymmetric case: " + rates + " rates " + str(dict_rates[rates]))
    plt.legend()
    fig.savefig('../figures/reachingtimes_asymm_' + rates + '.png')
    plt.close()  



def plot_holding(rates, zealots, times_z, contrarians, times_c):

    fig = plt.figure(figsize=(6,6))
    p1 = [2 * value for value in zealots]
    plt.plot(p1, np.log(times_z), 'b', linewidth = 1.5, label = 'Zealots')
    plt.plot(p1, np.log(times_c), 'r', linewidth = 1.5, label = 'Contrarians')
    plt.xlim(0,100)
    plt.ylim(-3, 10)
    plt.xlabel("Amount of disruptive individuals")
    plt.ylabel("Expected time to hold consensus in log scale")
    plt.title("Asymmetric case: " + rates + " rates " + str(dict_rates[rates]))
    plt.legend()
    fig.savefig('../figures/holdingtimes_asymm_' + rates + '.png')
    plt.close()  



# compute expected times for reaching consensus for all asymmetric rate pairs
# for rates in dict_rates:
#     verify_reaching_prism_zealots(rates)
#     verify_reaching_prism_contrarians(rates)
#     zealots, reachingtimes_z = read_reaching_prism_zealots(rates)
#     contrarians, reachingtimes_c = read_reaching_prism_contrarians(rates)
#     plot_reaching(rates, zealots, reachingtimes_z, contrarians, reachingtimes_c)

# compute expected times for holding consensus for all asymmetric rate pairs
for rates in dict_rates:
    verify_holding_prism_zealots(rates)
    verify_holding_prism_contrarians(rates)
    zealots, times_z = read_holding_prism_zealots(rates)
    contrarians, times_c = read_holding_prism_contrarians(rates)
    plot_holding(rates, zealots, times_z, contrarians, times_c)

