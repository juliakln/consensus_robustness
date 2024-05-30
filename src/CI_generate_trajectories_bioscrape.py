import bioscrape as bs
from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model

import numpy as np
import pylab as plt
import pandas as pd


# Read SBML Model
M = Model(sbml_filename = './models/CI_model_sbml2.xml')

timepoints = np.linspace(0, 250, 500)

num_trajectories = 250
exp_data = pd.DataFrame()
exp_data['timepoints'] = timepoints

for i in range(num_trajectories):
    result = py_simulate_model(timepoints, Model = M, stochastic = True)
    exp_data['X' + str(i)] = result['X']
    exp_data['Y' + str(i)] = result['Y']
    exp_data['U' + str(i)] = result['U']
    exp_data['Zx' + str(i)] = result['Zx']
    exp_data['Zy' + str(i)] = result['Zy']
    plt.plot(timepoints, exp_data['X' + str(i)], alpha=0.3)
    plt.plot(timepoints, exp_data['Y' + str(i)], alpha=0.3)
    plt.plot(timepoints, exp_data['U' + str(i)], alpha=0.3)
#plt.legend()
plt.xlabel('Time')
plt.ylabel('Species')
plt.show()

exp_data.to_csv('./data/CI2_trajectories_250.csv')
