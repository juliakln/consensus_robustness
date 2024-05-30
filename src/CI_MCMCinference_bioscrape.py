import matplotlib.pyplot as plt
import matplotlib as mpl
#%config InlineBackend.figure_f.ormats=['svg']
color_list = ['r', 'k', 'b','g','y','m','c']
mpl.rc('axes', prop_cycle=(mpl.cycler('color', color_list) ))
mpl.rc('xtick', labelsize=12) 
mpl.rc('ytick', labelsize=12)
import numpy as np
import bioscrape as bs
from bioscrape.types import Model
from bioscrape.inference import py_inference
from bioscrape.simulator import py_simulate_model
from bioscrape.types import Model
from bioscrape.sbmlutil import import_sbml
import pylab as plt
import pandas as pd

# Load model
#M_loaded_sbml = import_sbml('CI_model_sbml.xml')
M = Model(sbml_filename = './models/CI_model_sbml2.xml')

num_trajectories = 250

# Import data from CSV
# Import a CSV file for each experiment run
exp_data = []
for i in range(num_trajectories):
    df = pd.read_csv('./data/CI2_trajectories_250.csv', usecols = ['timepoints', 'X'+str(i), 'Y'+str(i), 'U'+str(i), 'Zx'+str(i), 'Zy'+str(i)])
    df.columns = ['timepoints', 'X', 'Y', 'U', 'Zx', 'Zy']
    exp_data.append(df)

prior = {'q1' : ['uniform', 0, 1]}

sampler, pid = py_inference(Model = M, exp_data = exp_data, measurements = ['X', 'Y', 'U', 'Zx', 'Zy'], time_column = ['timepoints'],
            nwalkers = 25, init_seed = "prior", nsteps = 1200, sim_type = 'stochastic',
            params_to_estimate = ['q1'], prior = prior)
report,_ = pid.plot_mcmc_results(sampler)

print('bla')


# use MCMC results to simulate model again and plot 1 trajectory
timepoints = pid.timepoints[0]
flat_samples = sampler.get_chain(discard=200, thin=15, flat=True)
inds = np.random.randint(len(flat_samples), size=1)
for ind in inds:
    sample = flat_samples[ind]
    for pi, pi_val in zip(pid.params_to_estimate, sample):
        M.set_parameter(pi, pi_val)
        res = py_simulate_model(timepoints, Model=M, stochastic = True)
    plt.plot(timepoints, res['X'], "C1", alpha=0.6, label='X')
    plt.plot(timepoints, res['Y'], "C2", alpha=0.6, label='Y')
    plt.plot(timepoints, res['U'], "C3", alpha=0.6, label='U')
plt.legend()

print('huwd')
# flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
# inds = np.random.randint(len(flat_samples), size=10)
# for ind in inds:
#     sample = flat_samples[ind]
#     plt.plot(exp_data['timepoints'], np.dot(np.vander(exp_data['timepoints'], 2), sample[:2]), "C1", alpha=0.1)
# plt.xlabel("x")
# plt.ylabel("y")