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
    ## analysis for stable consensus
#    analyse_stable_100('z')
    analyse_stable_500('z')
#    analyse_stable_1000('z')
    
#    analyse_stable_100('c')
    analyse_stable_500('c')
#    analyse_stable_1000('c')

#    analyse_stable_100_both()
#    analyse_stable_500_both()
#    analyse_stable_1000_both()

    ## analysis for switching consensus
#    analyse_switch_100('z')
    analyse_switch_500('z')
#    analyse_switch_1000('z')

#    analyse_switch_100('c')
    analyse_switch_500('c')
#    analyse_switch_1000('c')

#    analyse_switch_100_both()
#    analyse_switch_500_both()
#    analyse_switch_1000_both()


if __name__ == "__main__":
    sys.exit(main())