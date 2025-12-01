# Robustness of Consensus Reaching in Robotic Swarms

This repository contains the models, scripts, and results used to study the robustness of consensus in robotic swarms. The work investigates how robustly swarms can reach and maintain consensus in the presence of disruptive individuals, specifically zealots and contrarians, under stochastic dynamics.
The underlying decision-making mechanisms are based on the Voter Model and the Cross-Inhibition Model, represented as chemical reaction networks (CRNs). We analyse both symmetric scenarios (options X and Y have equal quality) and asymmetric scenarios, where one option (X) has superior quality. 
Properties are expressed in Bounded Linear Temporal Logic (BLTL) and evaluated using Plasmalab. Python scripts automate the execution of the Plasmalab command-line interface.


**Published results:**  

[Quantifying consensus in stochastic swarms with disruptive individuals](https://ieeexplore.ieee.org/abstract/document/11186823/)  
Julia Klein and Tatjana Petrov  
*European Control Conference (ECC)*, 272-277, 2025. 

[Exploring Consensus Robustness in Swarms with Disruptive Individuals](https://link.springer.com/chapter/10.1007/978-3-031-75107-3_3)  
Julia Klein, Alberto d’Onofrio, and Tatjana Petrov  
*International Symposium on Leveraging Applications of Formal Methods*, 33-48, 2024.

## Repository structure

### data/
Datasets for experiments and simulations  
./histograms/  
Distribution of x over different time steps for symmetric and asymmetric model

### figures/
Figures and plots generated from experiments

### inference_results/
Outcome files of model checking

### models/
PRISM and Plasmalab models and properties

### src/
Scripts and analysis tools in Python




# Running Plasmalab

- *cd 'plasmalab-1.4.5-SNAPSHOT/demos*
- *../plasmagui.sh launch*

# Expected times 
- use PRISM, in options: change "Engine" to 'sparse', increase value of "Termination max. iterations"; then verify property
- models:
  - zealots: models/tanja_cross_inhibition.prism 
  - contrarians: models/tanja_cross_inhibition_contrarians.prism 
  - both: models/tanja_cross_inhibition_z_and_c.prism
- properties:
  - models/tanja_CI.props 

### Reaching time
- reward: !(((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1; --> for not being in a majority state
- property: R=? [ F ("x_wins"|"y_wins") ]; --> for reaching a majority state

### Holding time
- reward: (((x+Zx>=m)&(x>=y+d)) | ((y+Zy>=m)&(y>=x+d))) : 1; --> for being in a majority state
- property:  R=? [ F (("x_wins"|"y_wins")&((!("y_wins")|(F (!("y_wins"))))&(!("x_wins")|(F (!("x_wins")))))) ]; --> for leaving a majority state after reaching it


