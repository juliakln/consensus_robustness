"""
Run all robustness analyses
- Probability of reaching a stable consensus over #stubborn individuals
    - Groups of 100, 500, and 1000 individuals
    - Model with zealots, contrarians, or both
- Probability of switching consensus over #stubborn individuals
    - Groups of 100, 500, and 1000 individuals
    - Model with zealots, contrarians, or both

Results are saved in figures/
"""

import os
import sys
import time

from analysis_stableconsensus import *
from analysis_switchingconsensus import *

# change to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# define file for model and property
model = "../models/consensus_model.rml"
property = "../models/consensus.bltl"



def main():

    # times = []

    analyse_stable(N=100, stubborn='c')
    analyse_switch(N=100, stubborn='c')

    # times.extend([analyse_stable(N=100, stubborn='z'), 
    #              analyse_stable(N=300, stubborn='z'),
    #              analyse_stable(N=500, stubborn='z'),
    #              analyse_stable(N=700, stubborn='z'),
    #              analyse_stable(N=1000, stubborn='z')])
    
    # print(times)

    # times.extend([analyse_stable(N=100, stubborn='c'),
    #               analyse_stable(N=300, stubborn='c'),
    #               analyse_stable(N=500, stubborn='c'),
    #               analyse_stable(N=700, stubborn='c'),
    #               analyse_stable(N=1000, stubborn='c')])

    # print(times)

    # times.extend([analyse_switch(N=100, stubborn='z'),  # 1413.57 seconds
    #               analyse_switch(N=300, stubborn='z'),  # 10436.64 seconds 
    #               analyse_switch(N=500, stubborn='z'),  # 75047.39 seconds
    #               analyse_switch(N=700, stubborn='z'),  # 95726.82 seconds 
    #               analyse_switch(N=1000, stubborn='z')])# 99357.25 seconds

    # print(times)

    # times.extend([analyse_switch(N=100, stubborn='c'),
    #               analyse_switch(N=300, stubborn='c'),
    #               analyse_switch(N=500, stubborn='c'),
    #               analyse_switch(N=700, stubborn='c'),
    #               analyse_switch(N=1000, stubborn='c')])

    # print(times)

if __name__ == "__main__":
    sys.exit(main())