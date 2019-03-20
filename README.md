# DiploSHIC-machine-learning

## simulation_withtimeserial.py:

uses to simulate population with population pickled at each timepoint from fwdpy11 v 0.2. Theta/rho in used is 11000 as we have hypothetical windows go from 0 to 11 in units while model simulations only go from 0-1 units. 

## final_sampling.py

Each population is then sampled with final_sampling.py and also converted to statistic fvec file by diploshic and to be used by diploshic. 

## concatenated_fvec.sh

this script concatenated all statistic from all 110 windows for each replicate and each time point. After this script is done, each resulted file is also has the header from original statistic fvec file append to the first line.

## predictionscript.py

This script is used to make predictions by diploSHIC for resulted files from concatenated_fvec.sh using models generated earlier in the project.

## getmeangeneticvalue.py

This script is used to record genetic value of each replicate for each timepoint. Then, all genetic value from the sample time point is averaged using R.

## Parallel_generation_of_model.sh

This is used to generate three prediction models with theta/rho = 1000 for each simulation with discoal. Alpha value is changing with 25_250/ 250_2500/2500_25000 using -Pa option. 
