# QPGP
Code for the manuscript: T-SP-30202-2023 - IEEE Transactions on Signal Processing -"Quasi-Periodic Gaussian Process Modeling of Pseudo-Periodic Signals"

------------------------------------------------------------------

Abstract: 
Abstract--- Pseudo-periodic signals are frequently encountered in modern scientific and engineering applications. Most current signal modeling methods focus on strictly periodic signals, and they fail to account for both the within- and between-period correlations of pseudo-periodic signals, which may lower the modeling and prediction accuracy. To address this challenge, we develop a novel quasi-periodic Gaussian process method for signals collected at grids. It can well model the within- and between-period correlations and has an easy-to-interpret structure that can quantify the degree of cycle oscillations in pseudo-periodic signals. Moreover, to speed up the model estimation, prediction, and simulation, we propose a fast composite likelihood approach that decomposes the full likelihood in an exact manner. The acceleration is achieved by leveraging the fast Fourier transform (FFT) and exploiting the Kronecker structure of the within- and between-period correlations. Its superior performance over state-of-the-art methods is demonstrated through a simulation study and two real case studies on sunspot period estimation and cardiac arrhythmia monitoring.

------------------------------------------------------------------

Instructions:
1. Please download and install MATLAB to use this code.
2. Add \routine\ as well as \smt\ package to your workpath.
3. fit_QPGP.m in the routine package is the main function to fit the QPGP model.

------------------------------------------------------------------

Notes:
1. The input data for fitting QPGP model is supposed to be grid.
2. routine package also contains the comparison methods CPGP, EPGP, FNLS, NRCPE and MLPE.
3. The results and figures in this paper can be reproduced by the code.



