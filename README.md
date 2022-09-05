# AquinoRocco_RepressedGene_ExtrinsicWhiteNoiseExample
This small code reproduces numerical results of the influence of extrinsic white noise on the degradation rate in the repressed gene system (Eqs 3.1-3.2, Figure 2) from the paper "Bimodality in gene expression without feedback: from Gaussian white noise to log-normal coloured noise" by Gerardo Aquino and Andrea Rocco (2020). 


Please note comments within the Matlab script (the code is not optimized). 

To run the C++ solver:

g++ -g -o SDEsolver Main.cpp SDE_methods.cpp -lm -std=c++11
./SDEsolver
