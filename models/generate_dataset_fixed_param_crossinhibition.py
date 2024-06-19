from cProfile import label
import sys
import warnings
import numpy as np
import stochpy
import pandas as pd
from tqdm import tqdm
import pickle
import time
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class AbstractionDataset(object):

    def __init__(self, n_init_states, n_trajs, state_space_bounds, model_name, time_step, T):
        # state_space_bounds : shape = (state_space_dim,2)
        # param_space_bounds : shape = (param_space_dim,2)
        self.n_init_states = n_init_states
        self.n_trajs = n_trajs
        self.n_training_points = n_init_states*n_trajs
        self.state_space_bounds = state_space_bounds
        self.state_space_dim = state_space_bounds.shape[0]
        self.stoch_mod = stochpy.SSA(IsInteractive=False)
        self.model_name = model_name
        self.directory_name = '/Users/juliaklein/Documents/robots/'
        self.stoch_mod.Model(model_name+'.psc', dir=self.directory_name)  # J: add directory
        self.time_step = time_step
        self.T = T # end time
        self.N = None # population size

    def set_popul_size(self, N):
        self.N = N


    def time_resampling(self, data):
        time_index = 0
        # new time array to transform to time grid (equally spaced)
        time_array = np.linspace(0, self.T, num=int(self.T / self.time_step+1))  # J: shape(17,)
        # new_data contains data with new time scan
        # 1. column contains instants of time -> corresponds to time_array
        new_data = np.zeros((time_array.shape[0], data.shape[1]))  
        new_data[:, 0] = time_array
        for j in range(len(time_array)):
            while time_index < data.shape[0] - 1 and data[time_index + 1][0] < time_array[j]:
                time_index = time_index + 1
            new_data[j, 1:] = data[time_index, 1:]      
        return new_data


    def set_initial_states(self, init_state):
        X = int(init_state[0])
        Y = int(init_state[1])
        U = int(init_state[2])
        #Zx = int(init_state[2])
        #Zy = int(init_state[3])
        self.stoch_mod.ChangeInitialSpeciesCopyNumber("X", X)
        self.stoch_mod.ChangeInitialSpeciesCopyNumber("Y", Y)
        self.stoch_mod.ChangeInitialSpeciesCopyNumber("U", U)
        #self.stoch_mod.ChangeInitialSpeciesCopyNumber("Zx", Zx)
        #self.stoch_mod.ChangeInitialSpeciesCopyNumber("Zy", Zy)

    def sample_initial_states(self, n_points=None):
        if n_points == None:
            n_points = self.n_init_states

        # J: could use state_space_bounds here
        #set_of_init_states = np.random.randint(low=30, high=200, size=(n_points, self.state_space_dim))
        #set_of_init_states = np.random.randint(low=100, high=101, size=(n_points, self.state_space_dim))
        #set_of_init_states = np.empty((n_points, self.state_space_dim))
        set_of_init_states = np.expand_dims(self.state_space_bounds[:,0], axis=0)  # NOTE: set initial values exactly to the values defined in .psc

        return set_of_init_states


    def generate_adjacency_matrix(self):
        # find out number of species and number of reactions
        num_species = self.state_space_dim
        num_reactions = len(self.stoch_mod.sim_rates_tracked)
        num_nodes = num_species + num_reactions

        A = np.zeros((num_nodes, num_nodes))
        # hard-coded stiochiometric coefficients => not possible to get them from PySCes model?
        # A[0,4] = 1
        # A[0,5] = 1
        # A[0,7] = 1
        # A[0,9] = 1
        # A[1,4] = 1
        # A[1,5] = 1
        # A[1,6] = 1
        # A[1,8] = 1
        # A[2,8] = 1
        # A[3,9] = 1
        # A[4,0] = 2
        # A[5,1] = 2
        # A[6,0] = 1
        # A[7,1] = 1
        # A[8,0] = 1
        # A[8,2] = 1
        # A[9,1] = 1
        # A[9,3] = 1
        with open("/Users/juliaklein/Documents/robots/"+self.model_name+"_adjacency_matrix.pickle", 'wb') as file:
            pickle.dump(A, file)


    def generate_training_set(self):

        # J: Y contains initial states
        Inits = np.zeros((self.n_training_points,self.state_space_dim))
        # J: X contains trajectories
        # J: I prefer the trajectories starting from time point 0 (initial state)
        # But we also save initial state in Y to have both options 
        #X = np.zeros((self.n_training_points, int(self.T/self.time_step), self.state_space_dim))
        Xinit = np.zeros((self.n_training_points, int(self.T/self.time_step)+1, self.state_space_dim))

        # J: sample different initial states 
        initial_states = self.sample_initial_states()

        # J: need count to store info for each trajectory and each initial state, 2 loops
        count, avg_time = 0, 0
        for i in tqdm(range(self.n_init_states)):
            self.set_initial_states(initial_states[i,:])
            for t in range(self.n_trajs):
                begin_time = time.time()
                #self.stoch_mod.DoStochSim(method="Direct", trajectories=self.n_trajs, mode="time", end=self.T)
                self.stoch_mod.DoStochSim(method="Direct", trajectories=1, mode="time", end=self.T)
                ntraj_time = time.time()-begin_time
                #print("Time for 1 traj: ", one_traj_time)

                # write SIR.psc_species_timeseries1.txt -> contains only current trajectory, later saved to pickle file
                self.stoch_mod.Export2File(analysis='timeseries', datatype='species', IsAverage=False, directory=self.directory_name, quiet=False)
                avg_time += ntraj_time

                datapoint = pd.read_table(filepath_or_buffer=self.directory_name+'/'+self.model_name+'.psc_species_timeseries1.txt', delim_whitespace=True, header=1).drop(labels="Reaction", axis=1).drop(labels='Fired', axis=1).drop("N",axis = 1).values

                new_datapoint = self.time_resampling(datapoint)
                # J: trajectories starting from time point 0 (initial state)
                Xinit[count,:,:] = new_datapoint[0:,1:self.state_space_dim+1]
                #X[count,:,:] = new_datapoint[1:,1:self.state_space_dim+1]
                Inits[count,:] = initial_states[i,:self.state_space_dim]

                count += 1
        print("average time to gen 1 traj with SSA: ", avg_time/self.n_trajs)
        self.Xinit = Xinit
        self.Y_s0 = Inits
        self.generate_adjacency_matrix()


    def generate_validation_set(self, n_init_states, n_trajs):

        Inits = np.zeros((n_init_states,self.state_space_dim))
        # J: trajectories starting from time point 0 (initial state) 
        #X = np.zeros((n_init_states,  n_trajs,int(self.T/self.time_step), self.state_space_dim))
        Xinit = np.zeros((n_init_states, n_trajs, int(self.T/self.time_step)+1, self.state_space_dim))

        initial_states = self.sample_initial_states()  
            
        for ind in range(n_init_states):
            # J: iterate through previously set number of initial states
            self.set_initial_states(initial_states[ind,:])

            Inits[ind,:] = initial_states[ind,:self.state_space_dim]  # J: save selected initial state in Ys 
                
            for k in range(n_trajs):  # J: iterate through number of trajectories for each initial state
                #if k%100 == 0:
                print((ind+1), "/", n_init_states, "------------------K iter: ", (k+1), "/", n_trajs)
                
                self.stoch_mod.DoStochSim(method="Direct", trajectories=1, mode="time", end=self.T)  # J: simulate 1 trajectory
                self.stoch_mod.Export2File(analysis='timeseries', datatype='species', IsAverage=False, directory=self.directory_name, quiet=False)

                datapoint = pd.read_table(filepath_or_buffer=self.directory_name+'/'+self.model_name+'.psc_species_timeseries1.txt', delim_whitespace=True, header=1).drop(labels="Reaction", axis=1).drop(labels='Fired', axis=1).drop("N",axis = 1).values
                # J: a datapoint = trajectory
                # J: here make it to time grid (because of continuous time)
                new_datapoint = self.time_resampling(datapoint)
                # J: trajectories starting from time point 0 (initial state) 
                #X[ind,k,:,:] = new_datapoint[1:,1:self.state_space_dim+1]
                Xinit[ind,k,:,:] = new_datapoint[0:,1:self.state_space_dim+1]
            
        self.Xinit = Xinit
        self.Y_s0 = Inits


    # J: save only number of S, I, and R in pickle file
    def save_dataset(self, filename):
        dataset_dict = {"X": self.Xinit, "Y_s0": self.Y_s0}
        with open(filename, 'wb') as handle:
            pickle.dump(dataset_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
 

def run_training(filename, n_trajs, n_init_states):

    time_step = 1
    n_steps = 250
    T = n_steps*time_step

    state_space_bounds = np.array([[40,40], [40,40], [0,0]]) #, [10,10], [10,10]])

    sir_dataset = AbstractionDataset(n_init_states, n_trajs, state_space_bounds, 'CrossInhibition', time_step, T)
    start_time = time.time()
    sir_dataset.generate_training_set()
    print("Time to generate the training set w fixed param =", time.time()-start_time)

    sir_dataset.save_dataset(filename)


def run_validation(filename, n_trajs, n_init_states):

    time_step = 0.5
    n_steps = 32    # J: more trajectories but shorter than in training set?
    T = n_steps*time_step

    # J: unnecessary?
    state_space_bounds = np.array([[30,200],[30,200]])

    sir_dataset = AbstractionDataset(n_init_states, n_trajs, state_space_bounds, 'SIR', time_step, T)
    start_time = time.time()
    sir_dataset.generate_validation_set(n_init_states, n_trajs)
    print("Time to generate the validation set w fixed param =", time.time()-start_time)

    sir_dataset.save_dataset(filename)



def plot_trajs(filename):
    # J: read files and plot trajectories
    with (open(filename, "rb")) as openfile:
        while True:
            try:
                objects = (pickle.load(openfile))
            except EOFError:
                break

    # J: plot trajectory data
    fig = plt.figure()
    if len(np.shape(objects['X']))==3:   # training: all trajs without splitting of initial states
        plt.plot((objects['X'][:,:,0].T + 10), label='Y', c='C0')
        plt.plot(objects['X'][:,:,1].T, label='U', c='C1')
        plt.plot((objects['X'][:,:,2].T + 10), label='X', c='C2') # NOTE: add #zealots --> need to find that automatically        #plt.plot(objects['X'][:,:,2].T, label='Zx', c='C2')
        #plt.plot(objects['X'][:,:,3].T, label='Zy', c='C3')
    elif len(np.shape(objects['X']))==4: # validation: trajs separated according to initial state
        for i in range(len(objects['X'])):
            plt.plot(objects['X'][i,:,:,0].T, label ='S', c='C0')
            plt.plot(objects['X'][i,:,:,1].T, label ='I', c='C1')
            plt.plot(objects['X'][i,:,:,2].T, label ='R', c='C2')
    plt.xlabel('t')
    plt.ylabel('X[t]')
    handles, labels = plt.gca().get_legend_handles_labels()
    _, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in np.sort(ids)]
    labels = [labels[i] for i in np.sort(ids)]
    plt.legend(handles, labels, loc='best')    
    plt.show(block=True)
    plt.savefig('/Users/juliaklein/Documents/robots/CrossInhibition_training_set_test_Z1.png')



def main():

    # J: run both functions
    run_training("/Users/juliaklein/Documents/robots/CrossInhibition_training_set_test_Z1.pickle", n_trajs = 1, n_init_states = 1)
    #run_validation("/Users/juliaklein/Documents/Inferring_CRNs/Data/SIR/SIR_validation_set_large.pickle", n_trajs = 50, n_init_states = 200)

    # J: plot both simulated datasets
    plot_trajs("/Users/juliaklein/Documents/robots/CrossInhibition_training_set_test_Z1.pickle")
    #plot_trajs("../Inferring_CRNs/Data/SIR/SIR_validation_set_large.pickle")


if __name__ == "__main__":
    sys.exit(main())