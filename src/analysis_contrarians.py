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

model = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/contrarian_model.rml"
property = "/Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/demos/Consensus/consensus.bltl"

N = 100

def write_property(contrarians, majority, distance, transient, holding):
    threshold = int((majority/100)*N)
    f = open(property, "w")
    f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x+Cx>="+str(threshold)+") & ((x+Cx)-(y+Cy)>="+str(distance)+")) | ((y+Cy>="+str(threshold) + ") & ((y+Cy)-(x+Cx)>="+str(distance)+"))))")
    f.close()

def write_property_switching(zealots, majority, distance, transient, holding):
    f = open(property, "w")
    f.write("F<="+str(transient)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    f.close()

def write_model(contrarians):
    X = int(1/2*(N-2*contrarians))
    f = open(model, "w")
    f.write("""ctmc

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

            [conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
            [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
            [conxc]    (x>0) & (Cx>0) & (Cy<N) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
            [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cy : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
            [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
            [conyc]    (y>0) & (Cy>0) & (Cx<N) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
            [conxx]    (Cx>2) & (Cy<(N-1)) -> Cx*Cx : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
            [conyy]    (Cy>2) & (Cx<(N-1)) -> Cy*Cy : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

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
            [conxa] true -> qy : true;
            [conxb] true -> qy : true;	
            [conxc] true -> qy : true;
            [conya] true -> qx : true;
            [conyb] true -> qx : true;
            [conyc] true -> qx : true;
            [conxx] true -> qy : true;
            [conyy] true -> qx : true;

        endmodule""")
    f.close()

def simulate_plasma(majority, distance, transient, holding):
    # vary number of zealots
    #contrarians = [0,2,4,6,8,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45] 
    #contrarians = np.linspace(1,35,35)
    contrarians = np.linspace(2,50,49)


    # simulate chain with Plasma and save paths
    for z in contrarians:
        write_model(int(z))
        #write_property(int(z), majority, distance, transient, holding)
        write_property_switching(int(z), majority, distance, transient, holding)
        result = "./plasmares_" + str(int(z)) + ".txt"
        pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"=4239 -f proba --progress -o " + result
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
    plt.xlabel('Amount of Contrarians Cx = Cy')
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
    dir_con = '/Users/juliaklein/Documents/consensus_robustness/inference_results/contrarians_switch_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    #simulate_plasma(majority, distance, transient, holding)
    #plot_probs(dir_con)
    probs = read_data(dir_con)

    return probs

def analyse_probabilites():
    baseline = consensus(majority = 50, distance = 10, transient = 35, holding = 40)
    majority_lower = consensus(majority = 35, distance = 10, transient = 35, holding = 40)
    majority_higher = consensus(majority = 65, distance = 10, transient = 35, holding = 40)
    distance_lower = consensus(majority = 50, distance = 1, transient = 35, holding = 40)
    distance_higher = consensus(majority = 50, distance = 20, transient = 35, holding = 40)
    transient_lower = consensus(majority = 50, distance = 10, transient = 20, holding = 40)
    transient_higher = consensus(majority = 50, distance = 10, transient = 50, holding = 40)
    holding_lower = consensus(majority = 50, distance = 10, transient = 35, holding = 25)
    holding_higher = consensus(majority = 50, distance = 10, transient = 35, holding = 55)

    fig = plt.figure(figsize=(6,6))
    plt.plot(majority_lower.keys(), majority_lower.values(), 'b--', label = 'm=35')
    plt.plot(majority_higher.keys(), majority_higher.values(), 'b:', label = 'm=65')
    plt.plot(distance_lower.keys(), distance_lower.values(), 'y--', label = 'd=1')
    plt.plot(distance_higher.keys(), distance_higher.values(), 'y:', label = 'd=20')
    plt.plot(transient_lower.keys(), transient_lower.values(), 'r--', label = 't=20')
    plt.plot(transient_higher.keys(), transient_higher.values(), 'r:', label = 't=50')
    plt.plot(holding_lower.keys(), holding_lower.values(), 'g--', label = 'h=25')
    plt.plot(holding_higher.keys(), holding_higher.values(), 'g:', label = 'h=55')
    plt.plot(baseline.keys(), baseline.values(), 'k', linewidth = 1.5, label = 'Baseline')
    plt.xlim(0,50)
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50], [0,10,20,30,40,50,60,70,80,90,100], rotation='horizontal')
    #plt.xlim(0,210)
    #plt.xticks([5,35,65,95,125,155,185], [10,70,130,190,250,310,370], rotation='horizontal')
    plt.title("Cross-Inhibition Model with Contrarians")    
    plt.xlabel('Amount of Contrarians C = Cx + Cy')
    plt.ylabel('Probability to reach consensus')
    plt.legend()
    fig.savefig('../all_probabilities_contrarians.png')
    plt.close()    


def analyse_switching():
    baseline = consensus(majority = 50, distance = 10, transient = 35, holding = 10)
    distance_lower1 = consensus(majority = 50, distance = 5, transient = 35, holding = 10)
    distance_lower2 = consensus(majority = 50, distance = 7, transient = 35, holding = 10)
    distance_higher1 = consensus(majority = 50, distance = 12, transient = 35, holding = 10)
    distance_higher2 = consensus(majority = 50, distance = 15, transient = 35, holding = 10)
    holding_lower1 = consensus(majority = 50, distance = 10, transient = 35, holding = 5)
    holding_lower2 = consensus(majority = 50, distance = 10, transient = 35, holding = 7)
    holding_higher1 = consensus(majority = 50, distance = 10, transient = 35, holding = 12)
    holding_higher2 = consensus(majority = 50, distance = 10, transient = 35, holding = 15)

    fig = plt.figure(figsize=(6,6))
    plt.plot(distance_lower1.keys(), distance_lower1.values(), 'y-', label = 'd=5')
    plt.plot(distance_lower2.keys(), distance_lower2.values(), 'y:', label = 'd=7')
    plt.plot(distance_higher1.keys(), distance_higher1.values(), 'y--', label = 'd=12')
    plt.plot(distance_higher2.keys(), distance_higher2.values(), 'y-.', label = 'd=15')
    plt.plot(holding_lower1.keys(), holding_lower1.values(), 'g-', label = 'h=5')
    plt.plot(holding_lower2.keys(), holding_lower2.values(), 'g:', label = 'h=7')
    plt.plot(holding_higher1.keys(), holding_higher1.values(), 'g--', label = 'h=12')
    plt.plot(holding_higher2.keys(), holding_higher2.values(), 'g-.', label = 'h=15')
    plt.plot(baseline.keys(), baseline.values(), 'k', linewidth = 1.5, label = 'Baseline')
    plt.xlim(0,50)
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50], [0,10,20,30,40,50,60,70,80,90,100], rotation='horizontal')
    plt.title("Cross-Inhibition Model with Contrarians")    
    plt.xlabel('Amount of contrarians C = Cx + Cy')
    plt.ylabel('Probability of consensus switching')
    plt.legend()
    fig.savefig('../all_switches_contrarians.png')
    plt.close()  


def main():
   #analyse_probabilites()
   analyse_switching()

if __name__ == "__main__":
    sys.exit(main())