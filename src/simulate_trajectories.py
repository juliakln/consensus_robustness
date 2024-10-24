""""
Run PRISM from terminal to simulate model and create trajectories 

"""

from collections import defaultdict
import sys
import os
os.chdir("/Applications/prism-4.7-src/prism/bin/")
import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np


def simulate_prism(num_trajs):
    # PRISM - ONLY RUN TO SIMULATE NEW DATA
    # uncertain parameter values for training data
    #p1 = np.linspace(0, 1, 15) 


    # generate paths: 10 seconds long, save as txt
    for i in range(0, num_trajs):
            os.system('./prism ../../../../Users/juliaklein/Documents/consensus_robustness/models/zealots.pm -simpath time=100,snapshot=0.05 \
                    ../../../../Users/juliaklein/Documents/consensus_robustness/data/trajectories/path_' + str(i) + '.txt -simpathlen 5000')


def read_trajectories():

    # List to store data from all files
    data_list = {}
    i=0

    # Loop through each file
    for file in glob.glob("../../../../Users/juliaklein/Documents/consensus_robustness/data/trajectories/*.txt"):
      # Read the file into a DataFrame using whitespace as the delimiter
        df = pd.read_csv(file, delim_whitespace=True)
        
        # Extract the relevant columns: x, y, u
        df_subset = df[['x', 'y', 'u']]
        
        # Append to the list
        data_list[i] = df_subset
        i+=1

    return data_list

def get_values(t, data_list):
    x_values = []
    y_values = []
    u_values = []

    # Loop through each trajectory (each DataFrame in the dictionary)
    for filename, df in data_list.items():
        # Get the first row's x value (index 0)
        x_values.append(df.iloc[t]['x'])
        y_values.append(df.iloc[t]['y'])
        u_values.append(df.iloc[t]['u'])

    return x_values, y_values, u_values

def plot_histograms(t, data_list):
    x_values, y_values, u_values = get_values(t, data_list)

    # Plot the histogram of x-values
    plt.figure(figsize=(8, 6))
    plt.hist(x_values, bins=10, alpha=0.75, color='red', edgecolor='black')

    time = str(int(t/20))

    # Add labels and title
    plt.xlim(-5, 70)
    plt.xlabel('x value')
    plt.ylabel('Frequency')
    plt.title('Distribution of x values at timepoint ' + time)
    plt.savefig('../../../../Users/juliaklein/Documents/consensus_robustness/figures/trajs_x_'+time+'.png')


    # Plot the histogram of y-values
    plt.figure(figsize=(8, 6))
    plt.hist(y_values, bins=10, alpha=0.75, color='blue', edgecolor='black')

    # Add labels and title
    plt.xlim(-5, 70)
    plt.xlabel('y value')
    plt.ylabel('Frequency')
    plt.title('Distribution of y values at timepoint ' + time)
    plt.savefig('../../../../Users/juliaklein/Documents/consensus_robustness/figures/trajs_y_'+time+'.png')

    
    # Plot the histogram of u-values
    plt.figure(figsize=(8, 6))
    plt.hist(u_values, bins=10, alpha=0.75, color='grey', edgecolor='black')

    # Add labels and title
    plt.xlim(-5, 40)
    plt.xlabel('u value')
    plt.ylabel('Frequency')
    plt.title('Distribution of u values at timepoint ' + time)
    plt.savefig('../../../../Users/juliaklein/Documents/consensus_robustness/figures/trajs_u_'+time+'.png')


# def simulate_bee_prism_pmc():
#     p1 = np.linspace(0, 1, 100) 
#     # simulate chain with Prism and save paths
#     for p in p1:
#         result = "/Users/juliaklein/Documents/uni/MasterThesis/data/dtmc_1_pmc/case_" + str(p) + ".txt"
#         resultfile = open(result ,"w")
#         resultfile.close()
#         resultfile = open(result, "r+")
#         prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/uni/MasterThesis/models/bee_3.pm /Users/juliaklein/Documents/uni/MasterThesis/models/bee_3_p.pctl -const p=" + str(p) + " -exportresults " + result
#         prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
#         resultfile.close()

# def simulate_bee_prism2(paths):
#     # PRISM - ONLY RUN TO SIMULATE NEW DATA
#     # uncertain parameter values for training data
#     X = np.zeros((100,2))
#     p1 = np.linspace(0.1, 1, 10) 
#     p2 = np.linspace(0.1, 1, 10)
#     X[:,0] = np.repeat(p1, 10)
#     X[:,1] = np.tile(p2, 10)
#     # simulate chain with Prism and save paths
#     for x in X:
#         p1 = x[0]
#         p2 = x[1]
#         result = "/Users/juliaklein/Documents/uni/MasterThesis/data/dtmc_2/case_" + str(round(p1,1)) + "_" + str(round(p2,2)) + ".txt"
#         resultfile = open(result ,"w")
#         resultfile.close()
#         resultfile = open(result, "r+")
#         prismcommand = "/Applications/prism-4.7-src/prism/bin/prism /Users/juliaklein/Documents/uni/MasterThesis/models/bee_3.pm /Users/juliaklein/Documents/uni/MasterThesis/models/bee_3_p.pctl -const p=" + str(p1) + ",q1=" + str(p2) + " -sim -simsamples " + str(paths) + " -exportresults " + result
#         prismprocess = subprocess.check_call(prismcommand, stdin=None, stdout=None , stderr=None, shell=True)
#         resultfile.close()


# def read_bee_prism():
#     # save number of satisfactions for each value of p and experiment
#     satisfactions = defaultdict(list)
#     # read outcomes for all parameter values and compute number of satisfactions
#     for dirpath, dirs, files in os.walk("/Users/juliaklein/Documents/uni/MasterThesis/data/dtmc_1"):
#         for file in files:
#             if file.startswith("case"):
#                 p = round(float((file.split("_")[1]).rsplit(".", 1)[0]), 4)
#                 with open(os.path.join(dirpath, file), 'r') as f:
#                     data = f.readline()
#                     for last_line in f:
#                         pass
#                     S = float(last_line)
#                     satisfactions[p].append(S)

#     # Training data
#     paramValueSet = []
#     paramValueOutputs = []
#     for key in sorted(satisfactions):
#         paramValueSet.append(key)
#         paramValueOutputs.append(satisfactions[key])

#     return np.array(paramValueSet).reshape(-1,1), (np.array(paramValueOutputs).reshape(-1,1))




def main():
    # simulate_prism(1000)
    # data = read_trajectories()
    
    # plot_histograms(0, data)
    # plot_histograms(200, data)
    # plot_histograms(400, data)
    # plot_histograms(600, data)
    # plot_histograms(800, data)
    # plot_histograms(1000, data)
    # plot_histograms(1200, data)
    # plot_histograms(1400, data)
    # plot_histograms(1600, data)
    # plot_histograms(1800, data)
    # plot_histograms(2000, data)



    ### plots for consensus paper - reaching times and holding times 

    data = {
    "amount": [2, 16, 30, 50, 70, 80, 82, 84, 86, 88, 90],
    "Zealots": [5.95, 7.28, 9.02, 10.57, 12.04, 27.82, 39.94, 64.95, 128.85, 374.04, 2975.68],
    "Contrarians": [6.07, 7.81, 6.89, 1.95, 0.63, 0.52, 0.51, 0.51, 0.49, 0.49, 0.46]
    }

    # data = {
    #     "amount": [12, 14, 16, 24, 26, 28, 30, 34, 50, 70, 90],
    #     "Zealots": [20686.51, 16368.28, 13047.85, 210.98, 47.71, 14.13, 5.46, 1.61, 0.48, 0.34, 0.22],
    #     "Contrarians": [283.57, 22.53, 4.03, 0.37, 0.31, 0.27, 0.25, 0.21, 0.16, 0.15, 0.14]
    # }

    plt.plot(data["amount"], np.log(data["Zealots"]), color='b', label='Zealots')
    plt.plot(data["amount"], np.log(data["Contrarians"]), color='r', label='Contrarians')
    plt.xlabel('Amount of disruptive individuals')  # Set x-axis label
    plt.ylabel('Time in sec., log scale')  # Set y-axis label for Zealots
    plt.legend()
    plt.show()

    # Create figure and axis objects
    fig, ax1 = plt.subplots()

    # Plot Zealots data on the first y-axis
    #ax1.plot(data["amount"], (data["Zealots"]), color='b', label='Zealots')
    ax1.plot(data["amount"], np.log(data["Zealots"]), color='b', label='Zealots')
    ax1.set_xlabel('Amount of disruptive individuals')  # Set x-axis label
    #ax1.set_ylabel('Time in sec.')  # Set y-axis label for Zealots
    ax1.set_ylabel('Time in sec., log scale')  # Set y-axis label for Zealots
    #ax1.set_ylim(0, 3500)
    #ax1.set_ylim(0, 23000)
    ax1.set_ylim(-3, 12)
    ax1.tick_params(axis='y', labelcolor='b')  # Set tick colors for Zealots y-axis

    # Create a second y-axis for Contrarians
    ax2 = ax1.twinx()
    #ax2.plot(data["amount"], (data["Contrarians"]), color='r', label='Contrarians')
    ax2.plot(data["amount"], np.log(data["Contrarians"]), color='r', label='Contrarians')
    #ax2.set_ylabel('Reaching time in sec.', color='r')  # Set y-axis label for Contrarians
    #ax2.set_ylim(0,9)
    #ax2.set_ylim(0,325)
    ax2.set_ylim(-3,12)
    ax2.tick_params(axis='y', labelcolor='r')  # Set tick colors for Contrarians y-axis

    ax1.legend(loc = "upper left")
    ax2.legend()
    # Show the plot
    plt.show()
    
    print("hehee")


if __name__ == "__main__":
    sys.exit(main())