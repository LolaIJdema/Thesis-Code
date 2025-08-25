# Thesis-Code

This repository includes all the R-code to supplement my thesis.

The internal functions file includes all the clustMD functions that need to be loaded into your environment to run my functions. 

Simulated data is the code I used to simulate data to test clustMD_S and clustMD_I.

ClustMD_S is the function for mixed-type variable model based clustering with a sample, and clustMD_I is the function for incremental mixed-type variable model based clustering. The clustMD_S function can be used when both continuous, ordinal and nominal variables are present in the data, whereas clustMD_I can only be used when there are only continuous and ordinal variables. 

The input arguments of both these functions are mostly the same as those of the clustMD function. For both clustMD_S and clustMD_I, SampSize determines the proportion of the observations that are used in the sample on which clustMD is fit (0.2 in all models in the thesis) and for clustMD_I the argument GroupProp determines the proportion of observations with the lowest dentities to by put into a seperate group in the incremental steps (0.01 in all models in the thesis). For clustMD_I, the number of Monte Carlo samples (MaxIter) does not have to be specified.

