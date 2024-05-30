import matplotlib.pyplot as plt
import matplotlib as mpl
#%config InlineBackend.figure_f.ormats=['svg']
color_list = ['r', 'k', 'b','g','y','m','c']
mpl.rc('axes', prop_cycle=(mpl.cycler('color', color_list) ))
mpl.rc('xtick', labelsize=12) 
mpl.rc('ytick', labelsize=12)
import numpy as np
from bioscrape.simulator import py_simulate_model
from bioscrape.types import Model

# Create a list of species names
species = ["X", "Y", "U", "Zx", "Zy"]

# Create a list of parameters
params = [("q1", .04), ("q2", .01), ("noise", 0.)]  # scaled by N!

# Create reaction tuples
rxn1 = (["X", "Y"], ["X", "U"], "massaction", {"k":"q1"})
rxn2 = (["X", "Y"], ["Y", "U"], "massaction", {"k":"q2"})
rxn3 = (["X", "U"], ["X", "X"], "massaction", {"k":"q1"})
rxn4 = (["Y", "U"], ["Y", "Y"], "massaction", {"k":"q2"})
rxn5 = (["Y"], ["U"], "massaction", {"k":"noise"})
rxn6 = (["U"], ["Y"], "massaction", {"k":"noise"})
rxn7 = (["X"], ["U"], "massaction", {"k":"noise"})
rxn8 = (["U"], ["X"], "massaction", {"k":"noise"})
rxn9 = (["Y"], ["X"], "massaction", {"k":"noise"})
rxn10 = (["X"], ["Y"], "massaction", {"k":"noise"})
rxn11 = (["Y", "Zx"], ["U", "Zx"], "massaction", {"k":"q1"})
rxn12 = (["X", "Zy"], ["U", "Zy"], "massaction", {"k":"q2"})
rxn13 = (["U", "Zx"], ["X", "Zx"], "massaction", {"k":"q1"})
rxn14 = (["U", "Zy"], ["Y", "Zy"], "massaction", {"k":"q2"})
rxns = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6, rxn7, rxn8, rxn9, rxn10, rxn11, rxn12, rxn13, rxn14]

# Create initial condition dictionary
x0 = {"X":40, "Y":40, "U":0, "Zx":10, "Zy":10}
x1 = {"X":35, "Y":35, "U":0, "Zx":15, "Zy":15}

# Instantiate Model object
M = Model(species = species, parameters = params ,reactions = rxns, initial_condition_dict = x0)

# Simulate Model deterministically
#timepoints = np.arange(0, 10, 0.01)
timepoints = np.arange(0, 60, 0.1)
#timepoints = np.arange(0, 360, 0.5)
#print(timepoints)
#results_det = py_simulate_model(timepoints, Model = M)

#Simulate the Model Stochastically
results_stoch = py_simulate_model(timepoints, Model = M, stochastic = True)

# Plot the results
plt.figure(figsize = (12, 4))
#plt.subplot(121)
plt.title("Cross Inhibition")
#plt.plot(timepoints, results_det["X"], label = "deterministic X")
plt.plot(timepoints, results_stoch["X"] + results_stoch["Zx"], label = "stochastic X")
#plt.plot(timepoints, results_det["Y"], label = "deterministic Y")
plt.plot(timepoints, results_stoch["Y"] + results_stoch["Zy"], label = "stochastic Y")
#plt.plot(timepoints, results_det["U"], label = "deterministic U")
plt.plot(timepoints, results_stoch["U"], label = "stochastic U")
plt.legend()
plt.xlabel("Time")
plt.show(block=True)


# Write Model to SBML
M.write_sbml_model('./models/CI_model_sbml2.xml')

#f = open('CI_model_sbml.xml')
#print("Bioscrape Model converted to SBML:\n", f.read())