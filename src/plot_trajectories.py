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
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

def plot_consensus():

    file_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z10_x.txt'
    file_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z10_x.txt'
    file_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z10_x.txt'
    file_x4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z30_x.txt'
    file_x5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z30_x.txt'
    file_x6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z30_x.txt'
    file_x7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z50_x.txt'
    file_x8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z50_x.txt'
    file_x9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z50_x.txt'

    file_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z10_u.txt'
    file_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z10_u.txt'
    file_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z10_u.txt'
    file_u4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z30_u.txt'
    file_u5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z30_u.txt'
    file_u6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z30_u.txt'
    file_u7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj1_1000_Z50_u.txt'
    file_u8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj2_1000_Z50_u.txt'
    file_u9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj3_1000_Z50_u.txt'

#    file_x = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesx.txt'
#    file_y = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesy.txt'
#    file_u = '/Users/juliaklein/Documents/consensus_robustness/inference_results/consensus_speciesu.txt'

    #zealots = 50 # per option

    files_x = [file_x1, file_x2, file_x3, file_x4, file_x5, file_x6, file_x7, file_x8, file_x9]
    files_u = [file_u1, file_u2, file_u3, file_u4, file_u5, file_u6, file_u7, file_u8, file_u9]
    zealots = [50, 50, 50, 150, 150, 150, 250, 250, 250] 
    lines = ['-', '--', ':', '-', '--', ':', '-', '--', ':']
    lines2 = ['-', '-', '-', '--', '--', '--', ':', ':', ':']
    cols = ['r', 'darkred', 'b', 'darkblue', 'g', 'darkgreen', 'y', 'gold', 'm', 'mediumorchid', 'k', 'grey', 
            'c', 'lightseagreen', 'lightcoral', 'forestgreen', 'navy', 'peru', 'orange']
    allx = []
    allu = []
    alltime = []

    for i in np.arange(len(files_x)):
        time = []
        x = []
        u = []        
        
        file1 = files_x[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            x.append(float(content.split(' ')[1]))
        fx.close()

        file2 = files_u[i]
        fu = open(file2, 'r')
        while True:
            content = fu.readline()
            if not content:
                break
            u.append(float(content.split(' ')[1]))
        fu.close()

#        x = [v + zealots[i] for v in x] # add amount of zealots
        # scale to 100
        x = [(v+zealots[i]) / 10 for v in x] # add amount of zealots
        u = [(v+zealots[i]) / 10 for v in u] # add amount of zealots

        alltime.append(time)
        allx.append(x)
        allu.append(u)


    fig = plt.figure(figsize=(12,10))
    for j in [2,5,8]: 

        # plt.plot(alltime[j], allx[j], 'r', linestyle = lines2[j], label = 'X'+str(zealots[j]))
        # plt.plot(alltime[j], allu[j], 'darkgrey', linestyle = lines2[j], label = 'U'+str(zealots[j]))
        plt.plot(alltime[j], allx[j], cols[j], linestyle = '-', label = 'X'+str(zealots[j]))
        plt.plot(alltime[j], allu[j], cols[j+1], linestyle = ':', label = 'U'+str(zealots[j]))
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))

    plt.ylim(0,100)
    plt.xlim(0,80)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_1000_3.png')
    plt.close()   


    fig = plt.figure(figsize=(12,10))
    for j in [0,1,2]: 

        plt.plot(alltime[j], allx[j], 'r', linestyle = lines[j], label = 'X')
        plt.plot(alltime[j], allu[j], 'darkgrey', linestyle = lines[j], label = 'U')
        plt.axhline(y=zealots[j], color = 'g', linestyle = lines[j], label = 'Zealots X')

    plt.ylim(0,100)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_1000.png')
    plt.close()   

#     time = []
#     x = []
# #    y = []
#     u = []

    # fx = open(file_x, 'r')
    # while True:
    #     content = fx.readline()
    #     if not content:
    #         break
    #     time.append(float(content.split(' ')[0]))
    #     x.append(float(content.split(' ')[1]))
    # fx.close()

    # fy = open(file_y, 'r')
    # while True:
    #     content = fy.readline()
    #     if not content:
    #         break
    #     y.append(float(content.split(' ')[1]))
    # fy.close()

    # fu = open(file_u, 'r')
    # while True:
    #     content = fu.readline()
    #     if not content:
    #         break
    #     u.append(float(content.split(' ')[1]))
    # fu.close()

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



def plot_compare10():

    file_100_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_1_x.txt'
    file_100_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_2_x.txt'
    file_100_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_3_x.txt'

    file_100_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_1_u.txt'
    file_100_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_2_u.txt'
    file_100_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z10_3_u.txt'

    file_1000_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_1_x.txt'
    file_1000_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_2_x.txt'
    file_1000_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_3_x.txt'

    file_1000_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_1_u.txt'
    file_1000_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_2_u.txt'
    file_1000_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_3_u.txt'

    files = [file_100_x1, file_100_x2, file_100_x3, file_100_u1, file_100_u2, file_100_u3, 
             file_1000_x1, file_1000_x2, file_1000_x3, file_1000_u1, file_1000_u2, file_1000_u3]

    zealots = [5,5,5,0,0,0,50,50,50,0,0,0] 
    scaling = [1,1,1,1,1,1,10,10,10,10,10,10]
    lines = ['-', '-', '-', ':', ':', ':', '-', '-', '-', ':', ':', ':']
 
    cols = ['deepskyblue', 'blue', 'darkblue', 'salmon', 'darkred', 'red', 
            'palegreen', 'green', 'darkgreen', 'lightyellow', 'yellow', 'gold']
    
    labels = ['X, N=100, #1', 'X, N=100, #2', 'X, N=100, #3', 'U, N=100, #1', 'U, N=100, #2', 'U, N=100, #3',
              'X, N=1000, #1', 'X, N=1000, #2', 'X, N=1000, #3', 'U, N=1000, #1', 'U, N=1000, #2', 'U, N=1000, #3']

    all_values = []
    all_time = []

    for i in np.arange(len(files)):
        time = []
        value = []
        
        file1 = files[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            value.append(float(content.split(' ')[1]))
        fx.close()

        # scale to 100
        value = [(v+zealots[i]) / scaling[i] for v in value] # add amount of zealots

        all_time.append(time)
        all_values.append(value)


    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files)): 
        plt.plot(all_time[i], all_values[i], cols[i], linestyle = lines[i], label = labels[i])
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))
    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('3 trajectories with 10% Zealots in the system')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_10zealots.png')
    plt.close()   



def plot_compare30():

    file_100_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_1_x.txt'
    file_100_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_2_x.txt'
    file_100_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_3_x.txt'

    file_100_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_1_u.txt'
    file_100_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_2_u.txt'
    file_100_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z30_3_u.txt'

    file_1000_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_1_x.txt'
    file_1000_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_2_x.txt'
    file_1000_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_3_x.txt'

    file_1000_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_1_u.txt'
    file_1000_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_2_u.txt'
    file_1000_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_3_u.txt'

    files = [file_100_x1, file_100_x2, file_100_x3, file_100_u1, file_100_u2, file_100_u3, 
             file_1000_x1, file_1000_x2, file_1000_x3, file_1000_u1, file_1000_u2, file_1000_u3]

    zealots = [15,15,15,0,0,0,150,150,150,0,0,0] 
    scaling = [1,1,1,1,1,1,10,10,10,10,10,10]
    lines = ['-', '-', '-', ':', ':', ':', '-', '-', '-', ':', ':', ':']
 
    cols = ['deepskyblue', 'blue', 'darkblue', 'salmon', 'darkred', 'red', 
            'palegreen', 'green', 'darkgreen', 'lightyellow', 'yellow', 'gold']
    
    labels = ['X, N=100, #1', 'X, N=100, #2', 'X, N=100, #3', 'U, N=100, #1', 'U, N=100, #2', 'U, N=100, #3',
              'X, N=1000, #1', 'X, N=1000, #2', 'X, N=1000, #3', 'U, N=1000, #1', 'U, N=1000, #2', 'U, N=1000, #3']

    all_values = []
    all_time = []

    for i in np.arange(len(files)):
        time = []
        value = []
        
        file1 = files[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            value.append(float(content.split(' ')[1]))
        fx.close()

        # scale to 100
        value = [(v+zealots[i]) / scaling[i] for v in value] # add amount of zealots

        all_time.append(time)
        all_values.append(value)


    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files)): 
        plt.plot(all_time[i], all_values[i], cols[i], linestyle = lines[i], label = labels[i])
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))
    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('3 trajectories with 30% Zealots in the system')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_30zealots.png')
    plt.close()   



def plot_compare50():

    file_100_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_1_x.txt'
    file_100_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_2_x.txt'
    file_100_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_3_x.txt'

    file_100_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_1_u.txt'
    file_100_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_2_u.txt'
    file_100_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_100_Z50_3_u.txt'

    file_1000_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_1_x.txt'
    file_1000_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_2_x.txt'
    file_1000_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_3_x.txt'

    file_1000_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_1_u.txt'
    file_1000_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_2_u.txt'
    file_1000_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_3_u.txt'

    files = [file_100_x1, file_100_x2, file_100_x3, file_100_u1, file_100_u2, file_100_u3, 
             file_1000_x1, file_1000_x2, file_1000_x3, file_1000_u1, file_1000_u2, file_1000_u3]

    zealots = [25,25,25,0,0,0,250,250,250,0,0,0] 
    scaling = [1,1,1,1,1,1,10,10,10,10,10,10]
    lines = ['-', '-', '-', ':', ':', ':', '-', '-', '-', ':', ':', ':']
 
    cols = ['deepskyblue', 'blue', 'darkblue', 'salmon', 'darkred', 'red', 
            'palegreen', 'green', 'darkgreen', 'lightyellow', 'yellow', 'gold']
    
    labels = ['X, N=100, #1', 'X, N=100, #2', 'X, N=100, #3', 'U, N=100, #1', 'U, N=100, #2', 'U, N=100, #3',
              'X, N=1000, #1', 'X, N=1000, #2', 'X, N=1000, #3', 'U, N=1000, #1', 'U, N=1000, #2', 'U, N=1000, #3']

    all_values = []
    all_time = []

    for i in np.arange(len(files)):
        time = []
        value = []
        
        file1 = files[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            value.append(float(content.split(' ')[1]))
        fx.close()

        # scale to 100
        value = [(v+zealots[i]) / scaling[i] for v in value] # add amount of zealots

        all_time.append(time)
        all_values.append(value)


    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files)): 
        plt.plot(all_time[i], all_values[i], cols[i], linestyle = lines[i], label = labels[i])
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))
    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('3 trajectories with 50% Zealots in the system')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_50zealots.png')
    plt.close()   


"""
    Plot trajectories (X and U) for N=1000 with different proportions of zealots
    Plot 1 traj for Z = 5%, 10%, 20%, 30%, 40%, 50%
"""
def plot_all_1000():

    file_1000_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_1_x.txt'
    file_1000_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_1_x.txt'
    file_1000_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z20_1_x.txt'
    file_1000_x4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_1_x.txt'
    file_1000_x5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z40_1_x.txt'
    file_1000_x6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_1_x.txt'

    file_1000_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_1_u.txt'
    file_1000_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z10_1_u.txt'
    file_1000_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z20_1_u.txt'
    file_1000_u4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z30_1_u.txt'
    file_1000_u5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z40_1_u.txt'
    file_1000_u6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z50_1_u.txt'

    files = [file_1000_x1, file_1000_x2, file_1000_x3, file_1000_x4, file_1000_x5, file_1000_x6, 
             file_1000_u1, file_1000_u2, file_1000_u3, file_1000_u4, file_1000_u5, file_1000_u6]

    zealots = [25,50,100,150,200,250,0,0,0,0,0,0] 
    scaling = 10
    lines = ['-', '-', '-', '-', '-', '-', ':', ':', ':', ':', ':', ':']
 
    cols = ['deepskyblue', 'salmon', 'palegreen', 'gold', 'plum', 'tan',
            'blue', 'red', 'green', 'yellow', 'darkviolet', 'darkgoldenrod']
    
    labels = ['X, 5%Z', 'X, 10%Z', 'X, 20%Z', 'X, 30%Z', 'X, 40%Z', 'X, 50%Z',
              'U, 5%Z', 'U, 10%Z', 'U, 20%Z', 'U, 30%Z', 'U, 40%Z', 'U, 50%Z',]

    all_values = []
    all_time = []

    for i in np.arange(len(files)):
        time = []
        value = []
        
        file1 = files[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            value.append(float(content.split(' ')[1]))
        fx.close()

        # scale to 100
        value = [(v+zealots[i]) / scaling for v in value] # add amount of zealots

        all_time.append(time)
        all_values.append(value)


    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files)): 
        plt.plot(all_time[i], all_values[i], cols[i], linestyle = lines[i], label = labels[i])
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))
    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('Trajectories for N=1000')
    plt.legend()
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_1000_comparison.png')
    plt.close()   

def rename_files(directory):
    # List all files in the directory
    files = os.listdir(directory)
    
    # Filter files that match the pattern "plasmares_X.txt" and sort them in descending order of X
    plasma_files = [file for file in files if file.startswith("plasmares_") and file.endswith(".txt")]
    plasma_files.sort(key=lambda x: int(x[10:-4]), reverse=True)
    
    for file in plasma_files:
        # Extract the current number X from the filename
        current_number = int(file[10:-4])
        
        # Calculate the new number Y
        new_number = 2 * current_number
        
        # Construct the new filename
        new_filename = f"plasmares_{new_number}.txt"
        
        # Construct the full file paths
        old_filepath = os.path.join(directory, file)
        new_filepath = os.path.join(directory, new_filename)
        
        # Rename the file
        os.rename(old_filepath, new_filepath)
        print(f"Renamed {file} to {new_filename}")




# Define the system of differential equations for the model with zealots
def system(t, y, q, z):
    x, u = y
    dx_dt = -q * x**2 + q * x * (u - z) + q * z * u
    du_dt = -2 * dx_dt
    return [dx_dt, du_dt]

"""
    Plot trajectories (X and U) for N=1000 with 5% of zealots
    Plot 10 trajectories, ODE solution and mean values
"""
def plot_all_1000_5():

    file_1000_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_1_x.txt'
    file_1000_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_2_x.txt'
    file_1000_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_3_x.txt'
    file_1000_x4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_4_x.txt'
    file_1000_x5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_5_x.txt'
    file_1000_x6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_6_x.txt'
    file_1000_x7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_7_x.txt'
    file_1000_x8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_8_x.txt'
    file_1000_x9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_9_x.txt'
    file_1000_x10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_10_x.txt'

    file_1000_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_1_u.txt'
    file_1000_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_2_u.txt'
    file_1000_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_3_u.txt'
    file_1000_u4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_4_u.txt'
    file_1000_u5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_5_u.txt'
    file_1000_u6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_6_u.txt'
    file_1000_u7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_7_u.txt'
    file_1000_u8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_8_u.txt'
    file_1000_u9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_9_u.txt'
    file_1000_u10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_1000_Z5_10_u.txt'

    files_x = [file_1000_x1, file_1000_x2, file_1000_x3, file_1000_x4, file_1000_x5, file_1000_x6,  file_1000_x7, file_1000_x8, file_1000_x9, file_1000_x10]
    files_u = [file_1000_u1, file_1000_u2, file_1000_u3, file_1000_u4, file_1000_u5, file_1000_u6,  file_1000_u7, file_1000_u8, file_1000_u9, file_1000_u10]

    zealots = 25
    scaling = 10


    allx = []
    allu = []
    alltime = []

    for i in np.arange(len(files_x)):
        time = []
        x = []
        u = []        
        
        file1 = files_x[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            x.append(float(content.split(' ')[1]))
        fx.close()

        file2 = files_u[i]
        fu = open(file2, 'r')
        while True:
            content = fu.readline()
            if not content:
                break
            u.append(float(content.split(' ')[1]))
        fu.close()

        # scale to 100
        x = [(v+zealots) / scaling for v in x] # add amount of zealots
        u = [v / scaling for v in u] # scale to 100

        alltime.append(time)
        allx.append(x)
        allu.append(u)


    # Define the common time grid
    common_time_grid = np.linspace(0, 70, 10001)

    # Interpolate each trajectory onto the common time grid
    interpolated_data_x = []
    interpolated_data_u = []
    for i in np.arange(len(allx)):
        interpolator = interp1d(alltime[i], allx[i], kind='linear', fill_value='extrapolate')
        interpolated_data_x.append(interpolator(common_time_grid))
        interpolator = interp1d(alltime[i], allu[i], kind='linear', fill_value='extrapolate')
        interpolated_data_u.append(interpolator(common_time_grid))

    # Calculate the mean trajectory
    mean_x = np.mean(interpolated_data_x, axis=0)
    mean_u = np.mean(interpolated_data_u, axis=0)

#    mean_x = np.mean(allx, axis=0)
#    mean_u = np.mean(allu, axis=0)

    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files_x)): 
        plt.plot(alltime[i], allx[i], 'b', linewidth = 0.5)
        plt.plot(alltime[i], allu[i], 'lightgrey', linewidth = 0.5)
#        plt.axhline(y=zealots[j], color = 'g', linestyle = lines2[j], label = 'Zealots X'+str(zealots[j]))
    plt.plot(common_time_grid, mean_x, 'darkblue', label = 'Mean X traj.')
    plt.plot(common_time_grid, mean_u, 'dimgrey', label = 'Mean U traj.')

    # Constants
    q = 1
    z = 2.5 # only of one kind
    # Initial conditions
    x0 = (100-2*z)/2 # this will also be equal to the y-decided individuals
    u0 = 0
    # Time span
    t_span = (0, 0.3)  # From t=0 to t=5
    t_eval = np.linspace(0, 0.3, 100)  # Evaluation points
    # Solve the system
    solution = solve_ivp(system, t_span, [x0, u0], args=(q, z), t_eval=t_eval)

    # Plot the solutions
    plt.plot(solution.y[0]+z, 'r', linestyle = '--', label='ODE x(t)')
    plt.plot(solution.y[1], 'orange', linestyle = '--', label='ODE u(t)')

    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('Trajectories for N=1000 with 5% zealots')
    plt.legend()
    plt.grid(True)
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_1000_5_overview.png')
    plt.close()   



"""
    Plot trajectories (X and U) for group size N with perc_Z% of zealots
    Plot 10 trajectories, ODE solution and mean values
"""
def plot_overview(N, perc_Z):

    file_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_1_x.txt'
    file_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_2_x.txt'
    file_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_3_x.txt'
    file_x4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_4_x.txt'
    file_x5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_5_x.txt'
    file_x6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_6_x.txt'
    file_x7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_7_x.txt'
    file_x8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_8_x.txt'
    file_x9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_9_x.txt'
    file_x10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_10_x.txt'

    file_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_1_u.txt'
    file_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_2_u.txt'
    file_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_3_u.txt'
    file_u4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_4_u.txt'
    file_u5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_5_u.txt'
    file_u6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_6_u.txt'
    file_u7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_7_u.txt'
    file_u8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_8_u.txt'
    file_u9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_9_u.txt'
    file_u10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_10_u.txt'

    files_x = [file_x1, file_x2, file_x3, file_x4, file_x5, file_x6,  file_x7, file_x8, file_x9, file_x10]
    files_u = [file_u1, file_u2, file_u3, file_u4, file_u5, file_u6,  file_u7, file_u8, file_u9, file_u10]

    # number of zealots per option
    zealots = (N * (perc_Z/100)) / 2
    scaling = N / 100


    allx = []
    allu = []
    alltime = []

    for i in np.arange(len(files_x)):
        time = []
        x = []
        u = []        
        
        file1 = files_x[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            x.append(float(content.split(' ')[1]))
        fx.close()

        file2 = files_u[i]
        fu = open(file2, 'r')
        while True:
            content = fu.readline()
            if not content:
                break
            u.append(float(content.split(' ')[1]))
        fu.close()

        # scale to 100
        x = [(v+zealots) / scaling for v in x] # add amount of zealots
        u = [v / scaling for v in u] # scale to 100

        alltime.append(time)
        allx.append(x)
        allu.append(u)


    # Define the common time grid
    common_time_grid = np.linspace(0, 70, 10001)

    # Interpolate each trajectory onto the common time grid
    interpolated_data_x = []
    interpolated_data_u = []
    for i in np.arange(len(allx)):
        interpolator = interp1d(alltime[i], allx[i], kind='linear', fill_value='extrapolate')
        interpolated_data_x.append(interpolator(common_time_grid))
        interpolator = interp1d(alltime[i], allu[i], kind='linear', fill_value='extrapolate')
        interpolated_data_u.append(interpolator(common_time_grid))

    # Calculate the mean trajectory
    mean_x = np.mean(interpolated_data_x, axis=0)
    mean_u = np.mean(interpolated_data_u, axis=0)

    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files_x)): 
        plt.plot(alltime[i], allx[i], 'b', linewidth = 0.5)
        plt.plot(alltime[i], allu[i], 'lightgrey', linewidth = 0.5)
    plt.plot(common_time_grid, mean_x, 'darkblue', label = 'Mean X traj.')
    plt.plot(common_time_grid, mean_u, 'dimgrey', label = 'Mean U traj.')

    # Constants
    q = 1
    z = perc_Z / 2 # only of one kind
    # Initial conditions
    x0 = (100-2*z)/2 # this will also be equal to the y-decided individuals
    u0 = 0
    # Time span
    t_span = (0, 0.7)  # From t=0 to t=5
    t_eval = np.linspace(0, 0.7, 100)  # Evaluation points
    # Solve the system
    solution = solve_ivp(system, t_span, [x0, u0], args=(q, z), t_eval=t_eval)

    # Plot the solutions
    plt.plot(solution.y[0]+z, 'r', linestyle = '--', label='ODE x(t)')
    plt.plot(solution.y[1], 'orange', linestyle = '--', label='ODE u(t)')

    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('Trajectories for N=' + str(N) + ' with ' + str(perc_Z) + '% zealots')
    plt.legend()
    plt.grid(True)
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_' + str(N) + '_' + str(perc_Z) + '_overview.png')
    plt.close()   


"""
    Plot trajectories (X and U) for group size N with perc_Z% of zealots
    Plot 10 trajectories, ODE solution and mean values
"""
def plot_overview_one(N, perc_Z):

    file_x1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_1_x.txt'
    file_x2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_2_x.txt'
    file_x3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_3_x.txt'
    file_x4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_4_x.txt'
    file_x5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_5_x.txt'
    file_x6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_6_x.txt'
    file_x7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_7_x.txt'
    file_x8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_8_x.txt'
    file_x9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_9_x.txt'
    file_x10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_10_x.txt'

    file_u1 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_1_u.txt'
    file_u2 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_2_u.txt'
    file_u3 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_3_u.txt'
    file_u4 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_4_u.txt'
    file_u5 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_5_u.txt'
    file_u6 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_6_u.txt'
    file_u7 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_7_u.txt'
    file_u8 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_8_u.txt'
    file_u9 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_9_u.txt'
    file_u10 = '/Users/juliaklein/Documents/consensus_robustness/data/traj_' + str(N) + '_Z' + str(perc_Z) + '_10_u.txt'

    files_x = [file_x1, file_x2, file_x3, file_x4, file_x5, file_x6,  file_x7, file_x8, file_x9, file_x10]
    files_u = [file_u1, file_u2, file_u3, file_u4, file_u5, file_u6,  file_u7, file_u8, file_u9, file_u10]

    # number of zealots per option
    zealots = (N * (perc_Z/100)) / 2
    scaling = N / 100


    allx = []
    allu = []
    alltime = []

    for i in np.arange(len(files_x)):
        time = []
        x = []
        u = []        
        
        file1 = files_x[i]
        fx = open(file1, 'r')
        while True:
            content = fx.readline()
            if not content:
                break
            time.append(float(content.split(' ')[0]))
            x.append(float(content.split(' ')[1]))
        fx.close()

        file2 = files_u[i]
        fu = open(file2, 'r')
        while True:
            content = fu.readline()
            if not content:
                break
            u.append(float(content.split(' ')[1]))
        fu.close()

        # scale to 100
        x = [(v+zealots) / scaling for v in x] # add amount of zealots
        u = [v / scaling for v in u] # scale to 100

        alltime.append(time)
        allx.append(x)
        allu.append(u)


    # Define the common time grid
    common_time_grid = np.linspace(0, 70, 10001)

    # Interpolate each trajectory onto the common time grid
    interpolated_data_x = []
    interpolated_data_u = []
    for i in np.arange(len(allx)):
        interpolator = interp1d(alltime[i], allx[i], kind='linear', fill_value='extrapolate')
        interpolated_data_x.append(interpolator(common_time_grid))
        interpolator = interp1d(alltime[i], allu[i], kind='linear', fill_value='extrapolate')
        interpolated_data_u.append(interpolator(common_time_grid))

    # Calculate the mean trajectory
    mean_x = np.mean(interpolated_data_x, axis=0)
    mean_u = np.mean(interpolated_data_u, axis=0)

    fig = plt.figure(figsize=(12,10))
    for i in np.arange(len(files_x)): 
        plt.plot(alltime[i], allx[i], 'b', linewidth = 0.5)
        plt.plot(alltime[i], allu[i], 'lightgrey', linewidth = 0.5)
    plt.plot(common_time_grid, mean_x, 'darkblue', label = 'Mean X traj.')
    plt.plot(common_time_grid, mean_u, 'dimgrey', label = 'Mean U traj.')

    # Constants
    q = 1
    z = perc_Z / 2 # only of one kind
    # Initial conditions
    x0 = (100-2*z)/2 # this will also be equal to the y-decided individuals
    u0 = 0
    # Time span
    t_span = (0, 0.7)  # From t=0 to t=5
    t_eval = np.linspace(0, 0.7, 100)  # Evaluation points
    # Solve the system
    solution = solve_ivp(system, t_span, [x0, u0], args=(q, z), t_eval=t_eval)

    # Plot the solutions
    plt.plot(solution.y[0]+z, 'r', linestyle = '--', label='ODE x(t)')
    plt.plot(solution.y[1], 'orange', linestyle = '--', label='ODE u(t)')

    plt.ylim(0,100)
    plt.xlim(0,70)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Number of individuals, scaled to 100')
    plt.title('Trajectories for N=' + str(N) + ' with ' + str(perc_Z) + '% zealots')
    plt.legend()
    plt.grid(True)
    fig.savefig('/Users/juliaklein/Documents/consensus_robustness/figures/trajectories_' + str(N) + '_' + str(perc_Z) + '_overview.png')
    plt.close()   




def main():
    #plot_consensus()
    #plot_switching()
    #plot_indecision()
    #plot_compare10()
    #plot_compare30()
    #plot_compare50()
    #plot_all_1000()

    # Replace 'your_directory_path' with the actual path to your folder
 #   directory_path = '/Users/juliaklein/Documents/consensus_robustness/inference_results/probs_switchconsensus2_z_1000_(50,150,35,10)'

    # Call the function
    #rename_files(directory_path)


    #plot_all_1000_5()

    plot_overview(N = 1000, perc_Z = 60)


if __name__ == "__main__":
    sys.exit(main())