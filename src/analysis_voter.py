"""
Invoke Plasmalab to simulate model and compute probability
using Monte Carlo
Plot results
"""

import os
import sys
import numpy as np
import subprocess
from collections import defaultdict
import matplotlib.pyplot as plt

model = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/voter_model.rml"
property = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/consensus.bltl"

def write_property(zealots, majority, distance, transient, holding):
    N = 100
    threshold = int((majority/100)*(N-2*zealots))
    f = open(property, "w")
    f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x>="+str(threshold)+") & (x-y>="+str(distance)+")) | ((y>="+str(threshold) + ") & (y-x>="+str(distance)+"))))")
    f.close()

def write_property_switching(zealots, majority, distance, transient, holding):
    N = 100
    threshold = int((majority/100)*(N-2*zealots))
    f = open(property, "w")
    f.write("F<="+str(transient)+" ((x-y>="+str(distance)+" & (true U<="+str(holding)+" (y-x>="+str(distance)+"))) | (y-x>="+str(distance)+" & (true U<="+str(holding)+" (x-y>="+str(distance)+"))))")
    f.close()

def write_model(zealots):
    N = 100
    X = int(1/2*(N-2*zealots))
    f = open(model, "w")
    f.write("""ctmc

        const int Zx = """ + str(zealots) + """;
        const int Zy = """ + str(zealots) + """;
        const int N = 100;

        module voter
            
            x : [0..N] init """ + str(X) + """;
            y : [0..N] init """ + str(X) + """;
            
            [cix] 	   (x>0) & (y>0) & (x<N) -> x*y : (x'=x+1) & (y'=y-1); // x+y -> x+x
            [ciy] 	   (x>0) & (y>0) & (y<N) -> x*y : (x'=x-1) & (y'=y+1); // x+y -> y+y
            [zeaxa]    (y>0) & (x<N)	 -> y*Zx : (y'=y-1) & (x'=x+1);	// y+Zx -> x+Zx
            [zeaya]    (x>0) & (y<N)	 -> x*Zy : (x'=x-1) & (y'=y+1);	// x+Zy -> y+Zy

        endmodule

        // base rates
        const double qx = 0.01; 
        const double qy = 0.01; 

        // module representing the base rates of reactions
        module base_rates
            
            [cix] true -> qx : true;
            [ciy] true -> qy : true;
            [zeaxa] true -> qx : true;	
            [zeaya] true -> qy : true;

        endmodule""")
    f.close()

def simulate_plasma(majority, distance, transient, holding):
    # vary number of zealots
    zealots = [0,2,4,6,8,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45] 

    # simulate chain with Plasma and save paths
    for z in zealots:
        write_model(z)
        write_property(z, majority, distance, transient, holding)
        #write_property_switching(z, majority, distance, transient, holding)
        result = "./plasmares_" + str(z) + ".txt"
        pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"=1060 -f proba --progress -o " + result
        prismprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

def plot_probs(dir_con):

    probs = read_data(dir_con)

    #p1 = [0,2,4,6,8,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45] 
    #probs = [1.,1.,1.,1.,0.99623,0.99151,0.98962,0.97736,0.94623,0.81792,0.64623,0.36321,0.12453,0.0217,0.00283,0.,0.,0.,0.,0.,0.]
    #probs = [1.,1.,1.,0.99811,0.99811,0.99245,0.9783,0.96226,0.92358,0.81132,0.63113,0.41321,0.22736,0.09906,0.03208,0.00849,0.,0.,0.,0.,0.]

 #   x, y = zip(*probs)
    fig = plt.figure(figsize=(6,6))
    plt.plot(probs.keys(), probs.values())
    plt.title(dir_con.split('/')[-1])
    plt.xlabel('Amount of Zealots Zx = Zy')
    plt.ylabel('Probability to reach consensus')
    fig.savefig('./probabilities.png')
    plt.close()

def read_data(dir = os.getcwd()):
    probabilities = {}

    for dirpath, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith("plasmares_"):
                zealots = int((file.split("_")[1]).split(".")[0])
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.read()
                    probabilities[zealots] = float(data.split('\n')[0])
                
    probs = dict(sorted(probabilities.items()))

    return probs

def consensus(majority, distance, transient, holding):
    dir_con = '/Users/juliaklein/Documents/consensus_robustness/inference_results/voter_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    simulate_plasma(majority, distance, transient, holding)
    plot_probs(dir_con)

def main():
    # consensus(majority = 50, distance = 10, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 7, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 5, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 3, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 13, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 15, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 20, transient = 45, holding = 10)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 7)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 5)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 3)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 13)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 15)
    consensus(majority = 50, distance = 10, transient = 35, holding = 40)
    consensus(majority = 50, distance = 1,  transient = 35, holding = 40)
    consensus(majority = 40, distance = 10, transient = 35, holding = 40)
    consensus(majority = 60, distance = 10, transient = 35, holding = 40)
    consensus(majority = 50, distance = 5,  transient = 35, holding = 40)
    consensus(majority = 50, distance = 15, transient = 35, holding = 40)
    consensus(majority = 50, distance = 10, transient = 25, holding = 40)
    consensus(majority = 50, distance = 10, transient = 45, holding = 40)
    consensus(majority = 50, distance = 10, transient = 35, holding = 30)
    consensus(majority = 50, distance = 10, transient = 35, holding = 60)

if __name__ == "__main__":
    sys.exit(main())