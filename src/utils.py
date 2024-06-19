"""
Helper functions 
- Write models in RML and properties in BLTL for Plasmalab
    - Cross-inhibition model with zealots, contrarians, or both
    - Property of reaching stable consensus or switching consensus
- Run Plasmalab and compute probability of property, save in txt file
- Read txt data
- Plot probabilities over number of stubborn individuals in the system for different settings
"""

import os
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict


# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model and property
model = "../../models/consensus_model.rml"
property = "../../models/consensus.bltl"

# define colours and linetypes for plotting
colours = ['b', 'b', 'b', 'y', 'y', 'r', 'r', 'g', 'g', 'm', 'm', 'c', 'c']
linetypes = ['-', '--', '--', ':', ':', '--', '--', ':', ':', '--', '--', ':', ':']


"""
Write cross-inhibition model with stubborn individuals for Plasmalab
    stubborn: z (zealots), c (contrarians)
    N: total population size
    init: initial number of stubborns for each opinion

    Output: RML file model
    NOTE: scale rates with N
"""
def write_model(stubborn, N, init):
    # compute remaining number of pure agents
    X = int(1/2 * (N - 2 * init))
    f = open(model, "w")
    if stubborn == 'z':
        f.write("""ctmc

            const int Zx = """ + str(init) + """;
            const int Zy = """ + str(init) + """;
            const int N = """ + str(N) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                u : [0..N] init 0;
                
                [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
                [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
                [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
                [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y
                [zeaxa]    (y>0) & (u<N)	 -> y*Zx : (y'=y-1) & (u'=u+1);		// y+Zx -> u+Zx
                [zeaxb]    (u>0) & (x<N)	 -> u*Zx : (x'=x+1) & (u'=u-1);		// u+Zx -> x+Zx
                [zeaya]    (x>0) & (u<N)	 -> x*Zy : (x'=x-1) & (u'=u+1);		// x+Zy -> u+Zy
                [zeayb]    (u>0) & (y<N)	 -> u*Zy : (y'=y+1) & (u'=u-1);		// u+Zy -> y+Zy

            endmodule

            // base rates
            const double qx = """ + str(1/N) + """; 
            const double qy = """ + str(1/N) + """; 

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

            endmodule""")
    elif stubborn == 'c':
        f.write("""ctmc

            const int N = """ + str(N) + """;
            const int cN = """ + str(2*init) + """;

            module cross_inhibition
                
                x : [0..N] init """ + str(X) + """;
                y : [0..N] init """ + str(X) + """;
                u : [0..N] init 0;
                Cx : [0..N] init """ + str(init) + """;
                Cy : [0..N] init """ + str(init) + """;
                
                [cix] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x) & (y'=y-1) & (u'=u+1); // x+y -> x+u
                [ciy] 	   (x>0) & (y>0) & (u<N) -> x*y : (x'=x-1) & (y'=y) & (u'=u+1); // x+y -> y+u
                [rx] 	   (x>0) & (x<N) & (u>0) -> x*u : (x'=x+1) & (u'=u-1);  	// u+x -> 2x 
                [ry]       (y>0) & (y<N) & (u>0) -> y*u : (y'=y+1) & (u'=u-1);		// u+y -> 2y

                [conxa]    (x>0) & (Cy>0) & (u<N) -> x*Cy : (x'=x-1) & (u'=u+1);		// x+Cy -> u+Cy
                [conxb]    (u>0) & (Cy>0) & (y<N) -> u*Cy : (u'=u-1) & (y'=y+1);		// u+Cy -> y+Cy
                [conxc]    (x>0) & (Cx>0) & (Cy<cN) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
                [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cx : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
                [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
                [conyc]    (y>0) & (Cy>0) & (Cx<cN) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
                [conxx]    (Cx>2) & (Cy<(cN-1)) -> (Cx*(Cx-1)/2) : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
                [conyy]    (Cy>2) & (Cx<(cN-1)) -> (Cy*(Cy-1)/2) : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

            endmodule

            // base rates
            const double qx = """ + str(1/N) + """; 
            const double qy = """ + str(1/N) + """; 

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
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()



"""
Write cross-inhibition model with zealots & contrarians for Plasmalab
    N: total population size
    zealots: initial number of zealots for each opinion
    contrarians: initial number of contrarians for each opinion

    Output: RML file model
    NOTE: scale rates with N
"""
def write_model_both(N, zealots, contrarians):
    X = int(1/2 * (N - (2 * contrarians) - (2 * zealots)))
    f = open(model, "w")
    f.write("""ctmc

        const int Zx = """ + str(zealots) + """;
        const int Zy = """ + str(zealots) + """;
        const int N = """ + str(N) + """;
        const int cN = """ + str(2*contrarians) + """;

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
            [conxc]    (x>0) & (Cx>0) & (Cy<cN) -> x*Cx : (Cx'=Cx-1) & (Cy'=Cy+1);		// x+Cx -> x+Cy
            [conya]    (y>0) & (Cx>0) & (u<N) -> y*Cx : (y'=y-1) & (u'=u+1);		// y+Cx -> u+Cx
            [conyb]    (u>0) & (Cx>0) & (x<N) -> u*Cx : (u'=u-1) & (x'=x+1);		// u+Cx -> x+Cx
            [conyc]    (y>0) & (Cy>0) & (Cx<cN) -> y*Cy : (Cy'=Cy-1) & (Cx'=Cx+1);		// y+Cy -> y+Cx
            [conxx]    (Cx>2) & (Cy<(cN-1)) -> Cx*(Cx-1)/2 : (Cx'=Cx-2) & (Cy'=Cy+2);		// Cx+Cx->Cy+Cy
            [conyy]    (Cy>2) & (Cx<(cN-1)) -> Cy*(Cy-1)/2 : (Cy'=Cy-2) & (Cx'=Cx+2);		// Cy+Cy->Cx+Cx

            [zeaconx]  (Cx>0) & (Cy<cN) -> Cx*Zx : (Cx'=Cx-1) & (Cy'=Cy+1);      // Cx+Zx -> Cy+Zx
            [zeacony]  (Cy>0) & (Cx<cN) -> Cy*Zy : (Cy'=Cy-1) & (Cx'=Cx+1);      // Cy+Zy -> Cx+Zy
            
        endmodule

        // base rates
            const double qx = """ + str(1/N) + """; 
            const double qy = """ + str(1/N) + """; 

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
            [conxx] true -> qy : true;
            *[conyy] true -> qx : true;
            [zeaconx] true -> qy : true;
            [zeacony] true -> qx : true;

        endmodule""")
    
    f.close()



"""
Write property for reaching stable consensus with stubborn individuals for Plasmalab
    N: total population size
    stubborn: z (zealots), c (contrarians) or b (both)
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes

    Output: BLTL file property
"""
def write_property_stableconsensus(N = 100, stubborn = 'z', majority = 50, distance = 10, transient = 35, holding = 40):
    # compute absolute number to reach majority
    threshold = int((majority / 100) * N)
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x+Zx>="+str(threshold)+") & (x-y>="+str(distance)+")) | ((y+Zy>="+str(threshold) + ") & (y-x>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x+Cx>="+str(threshold)+") & ((x+Cx)-(y+Cy)>="+str(distance)+")) | ((y+Cy>="+str(threshold) + ") & ((y+Cy)-(x+Cx)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(transient)+" (G<="+str(holding)+" (((x+Cx+Zx>="+str(threshold)+") & ((x+Cx+Zx)-(y+Cy+Zy)>="+str(distance)+")) | ((y+Cy+Zy>="+str(threshold) + ") & ((y+Cy+Zy)-(x+Cx+Zx)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()


"""
Write property for switching consensus with stubborn individuals for Plasmalab
    stubborn: z (zealots), c (contrarians) or b (both)
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes

    Output: BLTL file property
"""
def write_property_switching(stubborn = 'z', distance = 10, transient = 35, holding = 40):
    f = open(property, "w")
    if stubborn == 'z':
        f.write("F<="+str(transient)+" (((x+Zx)-(y+Zy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Zy)-(x+Zx)>="+str(distance)+"))) | ((y+Zy)-(x+Zx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Zx)-(y+Zy)>="+str(distance)+"))))")
    elif stubborn == 'c':
        f.write("F<="+str(transient)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    elif stubborn == 'b':
        f.write("F<="+str(transient)+" (((x+Cx)-(y+Cy)>="+str(distance)+" & (true U<="+str(holding)+" ((y+Cy)-(x+Cx)>="+str(distance)+"))) | ((y+Cy)-(x+Cx)>="+str(distance)+" & (true U<="+str(holding)+" ((x+Cx)-(y+Cy)>="+str(distance)+"))))")
    else:
        raise Exception('Type of stubborn individual not supported.')
    f.close()



"""
Run Plasmalab with commandline and compute probability of stable/switch property using Monte Carlo
    stubborn: z (zealots), c (contrarians)
    N: total population size
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes
    range: range of stubborn individuals to compute robustness
    filename: where to save the result files
    samples: how many samples to use for Monte Carlo

    Return: dictionary with keys = #stubborns and values = probabilities, read from result files
"""
def stableconsensus(stubborn, N, majority, distance, transient, holding, range, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_model(stubborn, N, int(s))
            write_property_stableconsensus(N, stubborn, majority, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs

def switchconsensus(stubborn, N, majority, distance, transient, holding, range, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for s in range:
        result = "./plasmares_" + str(int(s)) + ".txt"
        if not os.path.exists(result):
            write_model(stubborn, N, int(s))
            write_property_switching(stubborn, distance, transient, holding)
            pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
            pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data(os.getcwd())
    os.chdir('../')

    return probs




"""
Run Plasmalab with commandline and compute probability of stable/switch property with both zealots & contrarians using Monte Carlo
    N: total population size
    majority: more than m% of population commits to same decision
    distance: difference of at least d between majority and those favouring opposite decision
    transient: consensus is achieved within t minutes from start of dynamics
    holding: group maintains consensus for at least h minutes
    range_z: range of zealots to compute robustness
    range_c: range of contrarians to compute robustness
    filename: where to save the result files
    samples: how many samples to use for Monte Carlo

    Return: dictionary with keys = #stubborns and values = probabilities, read from result files
"""
def stableconsensus_both(N, majority, distance, transient, holding, range_z, range_c, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for z in range_z:
        for c in range_c:
            result = "./plasmares_" + str(int(z)) + "_" + str(int(c)) + ".txt"
            if not os.path.exists(result):
                write_model_both(N, int(z), int(c))
                write_property_stableconsensus(N, 'b', majority, distance, transient, holding)
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
                pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data_2dim(dir_con)

    return probs

def switchconsensus_both(N, majority, distance, transient, holding, range_z, range_c, filename, samples):
    dir_con = '../inference_results/' + filename + '_('+str(majority)+','+str(distance)+','+str(transient)+','+str(holding)+')'
    if not os.path.exists(dir_con):
        os.makedirs(dir_con)
    os.chdir(dir_con)

    for z in range_z:
        for c in range_c:
            result = "./plasmares_" + str(int(z)) + "_" + str(int(c)) + ".txt"
            if not os.path.exists(result):
                write_model_both(N, int(z), int(c))
                write_property_switching('b', distance, transient, holding)
                pcommand = "sh /Users/juliaklein/Documents/Sonstiges/plasmalab-1.4.5-SNAPSHOT/plasmacli.sh launch -m "+model+":rml -r "+property+":bltl  -a montecarlo -A\"Total samples\"="+str(samples)+" -f proba --progress -o " + result
                pprocess = subprocess.check_call(pcommand, stdin=None, stdout=None , stderr=None, shell=True)

    probs = read_data_2dim(dir_con)

    return probs



"""
Plot robustness analysis
    stubborn: z (zealots), c (contrarians)
    results: dictionaries of probabilities for different consensus settings
    labels: specification of settings for plotting labels
    zealots: list of zealots for which we computed the probabilities
    ylabel: What to write on y axis
    figname: how to name file

    Output: lineplot png of probabilities over stubborn individuals for different property settings
"""
def plot_results(stubborn, results, labels, zealots, ylabel, figname):
    fig = plt.figure(figsize=(6,6))
    for i in range(1, len(results)):
        plt.plot(results[i].keys(), results[i].values(), linestyle = linetypes[i], color = colours[i], label = labels[i])
    plt.plot(results[0].keys(), results[0].values(), 'k', linewidth = 1.5, label = 'Baseline')
    #plt.xlim(0, 2*int(zealots[-1]))
    # get every 3rd element for labels
    zealots_labels = zealots[3 - 1::3]
    plt.xticks(zealots_labels, np.multiply(zealots_labels, 2).astype(int), rotation='horizontal')
    if stubborn == 'z':
        plt.xlabel('Amount of Zealots Z = Zx + Zy')
    elif stubborn == 'c':
        plt.xlabel('Amount of Contrarians C = Cx + Cy')
    else:
        raise Exception('Type of stubborn individual not supported.')
    plt.ylabel(ylabel)
    plt.legend()
    fig.savefig('../figures/' + figname + '.png')
    plt.close()   


""" 
Plot 2dim robustness analysis for model with both zealots & contrarians
    results: dictionaries of probabilities for different consensus settings
    range_z: list of zealots for which we computed the probabilities
    range_c: list of zealots for which we computed the probabilities
    figname: how to name file

    Output: contput plot png of probabilities over both zealots & contrarians for 1 setting
"""
def plot_results_2dim(results, range_z, range_c, figname):
    x, y = np.meshgrid(range_z, range_c)
    p = np.array(list(results.values())).reshape(25,-1)

    fig = plt.figure(figsize=(6,6))
    plt.contourf(x, y, p.T, 50, vmin=0, vmax=1)
    plt.xticks([5,10,15,20,25,30,35,40,45], [10,20,30,40,50,60,70,80,90], rotation='horizontal')
    plt.yticks([5,10,15], [10,20,30], rotation='horizontal')
    plt.xlabel('Amount of zealots Z = Zx + Zy')
    plt.ylabel('Amount of contrarians C = Cx + Cy')
    plt.colorbar()
    plt.clim(0,1)
    fig.savefig('../figures/' + figname + '.png')
    plt.close()



"""
Read txt files containing probability of satisfying property
    dir: current directory of txt files

    Output: dictionary of probabilities, sorted by #stubborn individuals
"""
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