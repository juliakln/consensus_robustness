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

model = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/both_model.rml"
property = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/consensus.bltl"

N = 100

def write_property(majority, distance, transient, holding):
    threshold = int((majority/100)*N)
    f = open(property, "w")
    f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x+Cx+Zx>="+str(threshold)+") & ((x+Cx+Zx)-(y+Cy+Zy)>="+str(distance)+")) | ((y+Cy+Zy>="+str(threshold) + ") & ((y+Cy+Zy)-(x+Cx+Zx)>="+str(distance)+"))))")
    f.close()

def write_property_switching(majority, distance, transient, holding):
    f = open(property, "w")
    f.write("F<="+str(transient)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    f.close()

def write_model(zealots, contrarians):
    X = int(1/2*(N-(2*contrarians)-(2*zealots)))
    f = open(model, "w")
    f.write("""ctmc

        const int Zx = """ + str(zealots) + """;
        const int Zy = """ + str(zealots) + """;
        const int N = """ + str(N) + """;

        module cross_inhibition
            
            x : [0..N] init """ + str(X) + """;
            y : [0..N] init """ + str(X) + """;
            u : [0..N] init 0;
            Cx : [0..N] init """ + str(contrarians) + """;
	        Cy : [0..N] init """ + str(contrarians) + """;
            
            [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
            [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
            [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
            [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

            [zeaxa]    (y>0) & (u<N)	 -> y*Zx : (y'=y-1) & (u'=u+1);		// y+Zx -> u+Zx
            [zeaxb]    (u>0) & (x<N)	 -> u*Zx : (x'=x+1) & (u'=u-1);		// u+Zx -> x+Zx
            [zeaya]    (x>0) & (u<N)	 -> x*Zy : (x'=x-1) & (u'=u+1);		// x+Zy -> u+Zy
            [zeayb]    (u>0) & (y<N)	 -> u*Zy : (y'=y+1) & (u'=u-1);		// u+Zy -> y+Zy

            [conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
            [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
            [conxc]    (x>0) & (Cx>0) & (Cy<N) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
            [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cy : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
            [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
            [conyc]    (y>0) & (Cy>0) & (Cx<N) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
            [conxx]    (Cx>2) & (Cy<(N-1)) -> Cx*Cx : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
            [conyy]    (Cy>2) & (Cx<(N-1)) -> Cy*Cy : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

            [zeaconx]  (Cx>0) & (Cy<N) -> Cx*Zx : (Cx'=Cx-1) & (Cy'=Cy+1);      // Cx+Zx -> Cy+Zx
            [zeacony]  (Cy>0) & (Cx<N) -> Cy*Zy : (Cy'=Cy-1) & (Cx'=Cx+1);      // Cy+Zy -> Cx+Zy
            
        endmodule

        // base rates
        const double qx = 0.01; 
        const double qy = 0.01; 

        // module representing the base rates of reactions
        module base_rates
            
            [cix] true -> qx : true;
            [ciy] true -> qy : true;
            [rx] true -> qx : true;
            [ry] true -> qy : true;
            [zeaxa] true -> qx : true;	
            [zeaxb] true -> qx : true;
            [zeaya] true -> qy : true;
            [zeayb] true -> qy : true;
            [conxa] true -> qy : true;
            [conxb] true -> qy : true;	
            [conxc] true -> qy : true;
            [conya] true -> qx : true;
            [conyb] true -> qx : true;
            [conyc] true -> qx : true;
            [zeaconx] true -> qy : true;
            [zeacony] true -> qx : true;

        endmodule""")
    
    f.close()

def simulate_plasma(majority, distance, transient, holding):
    # vary number of zealots
    #zealots = [0,2,4,6,8,10,12,14,16,18,20] 
    #contrarians = [0,2,4,6,8,10,12,14,16,18,20] 
    # zealots = np.linspace(0,30,31)
    # contrarians = np.linspace(0,20,21)
    # zealots = np.linspace(1,20,13)
    # contrarians = np.linspace(1,25,13)
    # zealots = np.linspace(1,19,10)
    # contrarians = np.linspace(1,35,18)
    zealots = np.linspace(1,49,25)
    contrarians = np.linspace(1,19,10)

    # simulate chain with Plasma and save paths
    for z in zealots:
        for c in contrarians:
            result = "./plasmares_" + str(int(z)) + "_" + str(int(c)) + ".txt"
            if not os.path.exists(result):
                write_model(int(z), int(c))
                write_property(majority, distance, transient, holding)
                #write_property_switching(majority, distance, transient, holding)
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"=1060 -f proba --progress -o " + result
                prismprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

def plot_probs(dir_con):

    probs = read_data(dir_con)

 #   x, y = zip(*probs)
    fig = plt.figure(figsize=(6,6))
    plt.plot(probs.keys(), probs.values())
    plt.title(dir_con.split('/')[-1])
    plt.xlabel('Amount of Contrarians Cx = Cy')
    plt.ylabel('Probability to reach consensus')
    fig.savefig('./probabilities.png')
    plt.close()


def plot_probs_2dim(dir_con):

    probs = read_data_2dim(dir_con)
    # y, x = np.meshgrid(np.linspace(0, 28, 29), np.linspace(0, 20, 21))
    # p = np.array(list(probs.values())).reshape(29,-1)
    x, y = np.meshgrid(np.linspace(1, 49, 25), np.linspace(1, 19, 10))
    p = np.array(list(probs.values())).reshape(25,-1)

    fig = plt.figure(figsize=(6,6))
    plt.contourf(x, y, p.T, 50, vmin=0, vmax=1)
#    plt.title(dir_con.split('/')[-1])
    plt.xticks([5,10,15,20,25,30,35,40,45], [10,20,30,40,50,60,70,80,90], rotation='horizontal')
    plt.yticks([5,10,15], [10,20,30], rotation='horizontal')
    plt.xlabel('Amount of zealots Z = Zx + Zy')
    plt.ylabel('Amount of contrarians C = Cx + Cy')
    plt.colorbar()
    plt.clim(0,1)
    fig.savefig('../all_probabilities_both.png')
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

def read_data_2dim(dir = os.getcwd()):
    probabilities = defaultdict(list)

    for dirpath, dirs, files in os.walk(dir):
        for file in files:
            if file.startswith("plasmares_"):
                zealots = int((file.split("_")[1]))
                contrarians = int((file.split("_")[2]).split(".")[0])
                with open(os.path.join(dirpath, file), 'r') as f:
                    data = f.read()
                    probabilities[(zealots, contrarians)] = float(data.split('\n')[0])
                
    probs = dict(sorted(probabilities.items()))

    return probs

def consensus(majority, distance, transient, holding):
    dir_con = '/Users/juliaklein/Documents/consensus_robustness/inference_results/both_baseline_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    simulate_plasma(majority, distance, transient, holding)
    plot_probs_2dim(dir_con)

def main():
    #consensus(majority = 50, distance = 10, transient = 45, holding = 10)
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
    # consensus(majority = 50, distance = 10, transient = 45, holding = 20)
    # consensus(majority = 40, distance = 10, transient = 35, holding = 40)
    # consensus(majority = 60, distance = 10, transient = 35, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 25, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 35, holding = 30)
    # consensus(majority = 50, distance = 10, transient = 35, holding = 60)
    #consensus(majority = 50, distance = 10, transient = 35, holding = 40)
    # consensus(majority = 50, distance = 1,  transient = 35, holding = 40)
    # consensus(majority = 40, distance = 10, transient = 35, holding = 40)
    # consensus(majority = 60, distance = 10, transient = 35, holding = 40)
    # consensus(majority = 50, distance = 5,  transient = 35, holding = 40)
    # consensus(majority = 50, distance = 15, transient = 35, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 25, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 45, holding = 40)
    # consensus(majority = 50, distance = 10, transient = 35, holding = 30)
    # consensus(majority = 50, distance = 10, transient = 35, holding = 60)

    # baseline switching
    consensus(majority = 50, distance = 10, transient = 35, holding = 10)


if __name__ == "__main__":
    sys.exit(main())