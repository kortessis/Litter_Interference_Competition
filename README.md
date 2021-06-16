# Litter_Interference_Competition
This is a repository for code associated with the theory of interference competition by way litter production.

The code presented here can be organized in three categories.
1. Code used to analyze, simulate, and present results on the continuous-time model and its variants found in the supplementary materials, section S5. 
2. Code used to analyze, simulate, and present results on the discrete-time annual-perennial model and its variants in the main text.
3. A Matlab live script of supplementary material S6 which estimates parameters using data. 
4. A colormap to recreate the figures in the correct color (colormap info here: https://bids.github.io/colormap/).




***** 1. Code for continuous-time model and variants *****

1a. Litter_Interference_Growth_Rates.m
This script simulates population growth in the model where resources are fixed and litter is the only regulatory factor. The script is specifically designed to recreate figure S5 of the supplement, although can be modified for more species and different parameter values.

Dependencies: viridis.m and Litter_ODE.m




1b. Litter_ODE.m
This is the ODE to pass to an ODE solver in matlab for the model with n plant species and litter.

Inputs to this function are:
1. r = w*a*S - m, a row vector of species maximum growth rates when L = 0. It has the same number of entries as species.
2. beta, a row vector of species sensitivities to litter. It has the same number of entries as there are species.
3. v, a row vector of species litter conversion rates
4. m, a row vector of species biomass death rates
5. c, a scalar for litter decomposition rate.

Dependencies: none




1c. Litter_Resouces_ODE.m
This is the ODE to pass to the ODE solver in matlab for the model with 2 plant species, litter, and resources. The function has the following inputs:
1. w, a vector of resource conversion rates
2. a, a vector of maximum resource uptake rates
3. beta, a vector of species sensitivity to litter
4. v, a vector of species biomass conversion to litter rates
5. S, the resource supply point
6. r, the resource replenishment rate when R < S.
7. c, litter decomposition rate

Dependencies: none




1d. Litter_ODE_Sims.m
Script for plotting the dynamics of the full continuous-time model where both litter and resources are dynamic. Parameters are distinguished by environmental parameters (S, c, and r), species responses to litter and resources (L* and R*), and species effects (m, v, and a). 

Figure 4 is created from this script.

Dependencies: viridis.m and Litter_Resources_ODE.m




***** 2. Code for discrete-time model and variants *****

2a. APL_Sim_Tree.m.
A function that simulates the annual-perennial-litter model given by eqns 9-14 in the main text. Inputs are 
1. gen, the number of timesteps to run the model
2. init_cond, a 4-d row vector with the initial values of each state variable in the order [NA(t=0),L(t=0),NS(t=0),NP(t=0)].
3. a cell-based array of parameters for the model, with the following structure
	cell 1, a 3 element vector with survival probabilities [sA, sP, p1, p2]
	cell 2, a 3 element vector of seed yield values, [yA, yP, f]
	cell 3, a 2 element vector of seed germination values, [gA, gP]
	cell 4, a 2 element vector of establishment probabilities, [eA, eP]
	cell 5, a 3 element vector of litter decay and production rates [bA, bP, d, BT]
	cell 6, a 3 element vector of competition coefficient, [alphaA, alphaP, gamma]
	cell 7, a 2 element vector of litter sensitivities, [betaA, betaP]

Output from this function is the a gen x 4 matrix with the values of the 4 dynamical variables evaluated for each integer time values.

Dependencies: none



2b. Coexistence_Analysis.m
A script that uses the analytical results to plot competitive outcomes in the phase space of lambdaP vs. lambdaA. This script reproduces figure 3 of the main text.

Dependencies: LitterEq.m and viridis.m.



2c. Community_Sim.m
A script to simulate dynamics in the multi-perennial version of the annual-perennial model. This script, as written, reproduces figures 7 (from the main text) and S2-S4 (from the supplementary materials). 

Dependencies: viridis.m, a color palette function.



2d. Model_Dynamics.m

Code that produces invasion boundaries in the phase space of annual and perennial maximum fitness. This script creates figure 2 in the main text. It is comprised of two parts. The first part shows how the coexistence region changes with the sign and magnitude of a tradeoff between species in sensitivity to litter, given fixed litter production. The second shows the effect of decomposition on the coexistence region. 

Dependencies: viridis.m and APL_BothLitter_Producers_Sim.m



2e. LitterEq.m

A function that calculates the single-species equilbrium litter density, eqns (S12) and (S24) in the supplement for the annual and perennial, respectively. Inputs are lambda_j, alpha_j, beta_j, the ratio of bj/d, and the ratio BT/d. For the perennail, \alpha is taken to be \alpha' = \alpha*(1 + \gamma*(1-p2)/p1). 

Dependencies: none



***** 3. Code for fitting the model to data *****

This is a Matlab live script (.mlx) giving the entire analysis. The identical code is presented as supplementary material S6.

Dependencies: viridis.m and LitterEq.m


***** 4. Ancillary code *****

4. viridis.m
This script is simply a color palette. Input an integer value and it provides matrix with three columns for RBG-color values to be used across the viridis color palette. 
***NOTE: I did not produce this colormap. The original colormap can be found here: https://bids.github.io/colormap/. 
