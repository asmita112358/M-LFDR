# M-LFDR
Code for MLFDR

1. EM_general_fun.R contains EM algorithm functions for normal alternative, i.e 4 component mixture model.
2. sim.R contains codes for simulations.
3. twostep_EM.R contains codes for two step EM, applicable for mixture normal alternative, recommended when the alternative has a high number of atoms.
4. Real1_pp.R is the code for preprocessing the first real data.
5. Real2_pp.R preprocesses the second real data, with survival outcome.
6. cox_fit.R fits a high dimensional cox regression model to the second real data, obtained from Real2_pp.R
