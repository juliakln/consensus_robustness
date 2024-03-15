import numpy as np
import gillespy2
import matplotlib.pyplot as plt
import pickle

class VoterModel(gillespy2.Model):
    def __init__(self, N, timesteps, endtime):
        # First call the gillespy2.Model initializer.
        super().__init__(self)

        # Define parameters for the rates of creation and dissociation.
        q1 = gillespy2.Parameter(name='q1', expression=1)
        q2 = gillespy2.Parameter(name='q2', expression=1)
        noise = gillespy2.Parameter(name='noise', expression=0.1)
        self.add_parameter([q1, q2, noise])
        self.N = N
        # Define variables for the molecular species representing M & D.
        X = gillespy2.Species(name='X', initial_value=100)
        Y = gillespy2.Species(name='Y', initial_value=100)
        #Zx = gillespy2.Species(name='Zx', initial_value=0)
        #Zy = gillespy2.Species(name='Zy', initial_value=0)
        #self.species_list = [X, Y, Zx, Zy]
        self.species_list = [X, Y]
        self.add_species(self.species_list)

        # The list of reactants and products for a Reaction object are
        # each a Python dictionary in which the dictionary keys are
        # Species objects and the values are stoichiometries of the
        # species in the reaction.
        r1 = gillespy2.Reaction(name="r1", rate=q1,
                                 reactants={X:1}, products={X:2})
        r2 = gillespy2.Reaction(name="r2", rate=q2,
                                 reactants={Y:1}, products={Y:2})
        r3 = gillespy2.Reaction(name="r3", rate=noise,
                         reactants={Y:1}, products={X:1})
        r4 = gillespy2.Reaction(name="r4", rate=noise,
                         reactants={X:1}, products={Y:1})
        self.add_reaction([r1, r2, r3, r4])

        # Set the timespan for the simulation.
        self.timespan(np.linspace(0, endtime, timesteps))


    def set_initial_state(self):
        #x = np.random.randint(self.N/2)
        #y = np.random.randint(self.N-x)
        #c = self.N-a-b
        x = 100
        y = 100
        self.species_list[0].initial_value = x
        self.species_list[1].initial_value = y
        #self.species_list[2].initial_value = c

N = 200
H = 15
model = VoterModel(N, timesteps=H, endtime=10)
npoints = 1#000
ntrajs = 100

trajectories = np.empty((npoints*ntrajs, H, 2))
c=0
for i in range(npoints):
    print(f'i={i+1}/{npoints}')
    model.set_initial_state()
    trajs = model.run(algorithm = "SSA",number_of_trajectories=ntrajs)
    for j in range(ntrajs):
        trajectories[c,:,0] = trajs[j]['X']
        trajectories[c,:,1] = trajs[j]['Y']
        #trajectories[c,:,2] = trajs[j]['C']
        c += 1

filename = f'voter_set_{npoints}x{ntrajs}.pickle'
data = {'trajs': trajectories}
with open(filename, 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)        
'''

fig = plt.figure()
for index in range(0, 5):
    trajectory = trajs[index]
    plt.plot(trajectory['time'], trajectory['A'], 'r')
    plt.plot(trajectory['time'], trajectory['B'],   'b')
    plt.plot(trajectory['time'], trajectory['C'],   'g')
plt.savefig('trajs/toy_test_-1.png')
'''