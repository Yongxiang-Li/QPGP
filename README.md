# QPGP
Code for the paper: Yongxiang LI, Yunji ZHANG, Qian XIAO, Jianguo WU*. (2023). Quasi-Periodic Gaussian Process Modeling of Pseudo-Periodic Signals. IEEE Transactions on Signal Processing. Accepted. DOI: 10.1109/TSP.2023.3316589
Link: https://ieeexplore.ieee.org/document/10256149?source=authoralert
------------------------------------------------------------------

Abstract: 
Pseudo-periodic signals are frequently encountered in modern scientific and engineering applications. Most current signal modeling methods focus on strictly periodic signals, and they may fail to account for both the within- and between-period correlations of pseudo-periodic signals, which could lower the modeling and prediction accuracy. To address this issue, we develop a novel quasi-periodic Gaussian process method for signals collected at grids. It can well model the within- and between-period correlations and has an easy-to-interpret structure that can quantify the magnitude of cycle oscillations in pseudo-periodic signals. To speed up the model estimation, prediction, and simulation, we further propose a fast composite likelihood approach that decomposes the full likelihood in an exact manner. This acceleration is achieved by leveraging the fast Fourier transform (FFT) and exploiting the Kronecker structure of the within- and between-period correlations. Its superior performance over some state-of-the-art methods is demonstrated through a simulation study and two real case studies on sunspot period estimation and cardiac arrhythmia monitoring.

------------------------------------------------------------------

Instructions:
1. Please download and install MATLAB to use this code.
2. Add \routine\ package to your workpath.
3. fit_QPGP.m in the routine package is the main function to fit the QPGP model.

------------------------------------------------------------------

Notes:
1. The input data for fitting QPGP model is supposed to be grid.
2. routine package also contains the comparison methods CPGP, EPGP, FNLS, NRCPE and MLPE.
3. The results and figures in this paper can be reproduced by the code.

------------------------------------------------------------------

Citation:
@article{Li2023QPGP,
  title={Quasi-Periodic Gaussian Process Modeling of Pseudo-periodic Signals},
  author={Li, Yongxiang and Zhang, Yunji and Xiao, Qian and Wu, Jianguo},
  journal={IEEE Transactions on Signal Processing},
  volume={71},
  number={12},
  pages={2700--2710},
  year={2023},
  publisher={IEEE}
}



