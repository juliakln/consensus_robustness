"""
Plot example trajectories
1) stable consensus
2) switching
3) persistent indecision 
"""

import os
import sys
import numpy as np
import subprocess
from collections import defaultdict
import matplotlib.pyplot as plt


def plot_consensus():

    file_x = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesx.txt'
    file_y = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesy.txt'
    file_u = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesu.txt'

    zealots = 15

    time = []
    x = []
    y = []
    u = []

    fx = open(file_x, 'r')
    while True:
        content = fx.readline()
        if not content:
            break
        time.append(float(content.split(' ')[0]))
        x.append(float(content.split(' ')[1]))
    fx.close()

    fy = open(file_y, 'r')
    while True:
        content = fy.readline()
        if not content:
            break
        y.append(float(content.split(' ')[1]))
    fy.close()

    fu = open(file_u, 'r')
    while True:
        content = fu.readline()
        if not content:
            break
        u.append(float(content.split(' ')[1]))
    fu.close()

    # add number of zealots to x and y for plotting
    allx = [vx+zealots for vx in x]
    ally = [vy+zealots for vy in y]

    fig = plt.figure(figsize=(6,6))
    plt.plot(time, allx, 'r', label = 'X + Zx')
    plt.plot(time, ally, 'b', label = 'Y + Zy')
    plt.plot(time, u, color = 'darkgrey', label = 'U')
    plt.ylim(0,100)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/inference_results/trajectories_consensus.png')
    plt.close()    

def plot_switching():

    file_x = '/Users/juliaklein/Documents/consensus_robustness/inference_results/switching_speciesx.txt'
    file_y = '/Users/juliaklein/Documents/consensus_robustness/inference_results/switching_speciesy.txt'
    file_u = '/Users/juliaklein/Documents/consensus_robustness/inference_results/switching_speciesu.txt'

    zealots = 15

    time = []
    x = []
    y = []
    u = []

    fx = open(file_x, 'r')
    while True:
        content = fx.readline()
        if not content:
            break
        time.append(float(content.split(' ')[0]))
        x.append(float(content.split(' ')[1]))
    fx.close()

    fy = open(file_y, 'r')
    while True:
        content = fy.readline()
        if not content:
            break
        y.append(float(content.split(' ')[1]))
    fy.close()

    fu = open(file_u, 'r')
    while True:
        content = fu.readline()
        if not content:
            break
        u.append(float(content.split(' ')[1]))
    fu.close()

    # add number of zealots to x and y for plotting
    allx = [vx+zealots for vx in x]
    ally = [vy+zealots for vy in y]

    fig = plt.figure(figsize=(6,6))
    plt.plot(time, allx, 'r', label = 'X + Zx')
    plt.plot(time, ally, 'b', label = 'Y + Zy')
    plt.plot(time, u, color = 'darkgrey', label = 'U')
    plt.ylim(0,100)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/inference_results/trajectories_switching.png')
    plt.close()    


def plot_indecision():

    file_x = '/Users/juliaklein/Documents/consensus_robustness/inference_results/indecision_speciesx.txt'
    file_y = '/Users/juliaklein/Documents/consensus_robustness/inference_results/indecision_speciesy.txt'
    file_u = '/Users/juliaklein/Documents/consensus_robustness/inference_results/indecision_speciesu.txt'

    zealots = 15

    time = []
    x = []
    y = []
    u = []

    fx = open(file_x, 'r')
    while True:
        content = fx.readline()
        if not content:
            break
        time.append(float(content.split(' ')[0]))
        x.append(float(content.split(' ')[1]))
    fx.close()

    fy = open(file_y, 'r')
    while True:
        content = fy.readline()
        if not content:
            break
        y.append(float(content.split(' ')[1]))
    fy.close()

    fu = open(file_u, 'r')
    while True:
        content = fu.readline()
        if not content:
            break
        u.append(float(content.split(' ')[1]))
    fu.close()

    # add number of zealots to x and y for plotting
    allx = [vx+zealots for vx in x]
    ally = [vy+zealots for vy in y]

    fig = plt.figure(figsize=(6,6))
    plt.plot(time, allx, 'r', label = 'X + Zx')
    plt.plot(time, ally, 'b', label = 'Y + Zy')
    plt.plot(time, u, color = 'darkgrey', label = 'U')
    plt.ylim(0,100)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/inference_results/trajectories_indecision.png')
    plt.close()    


def main():
    #plot_consensus()
    #plot_switching()
    plot_indecision()

if __name__ == "__main__":
    sys.exit(main())