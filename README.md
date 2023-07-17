# Formulation-of-CC-DCOPF
This set of codes implements the big-M theory, the improved big-M theory, and the GMM-based theory for the CC-DCOPF problem.

# Introduction
The code of M1, according to Ref. [1], is implemented in M1.m.
The code of M2, according to Ref. [2], is implemented in M2.m.
The code of M3, according to Ref. [3], is implemented in M3_1.m and M3_2.m.
For M3, the uncertainty margins are obtained in M3_1.m, which is loaded in M3_2.m for CC-DCOPF calculation.
case118WashingtonData.m and case118.m provide the system parameter of the IEEE 118-bus system.
nearoptimalsolutions.mat is the near-optimal solution for 1000 scenarios obtained by the proposed method.
case118_uncertaintydata_100.mat and case118_uncertaintydata_1000.mat are the uncertainty data of the IEEE 118-bus system for 100 scenarios and 1000 scenarios, respectively.
mixGaussEm.m is for GMM theory.

# Matlab Codes Environment
Recommand Matlab Version: MATLAB R2020a
Required packages: Yalmip, Matpower
Required solver: Gurobi

# References
[1]	Qiu, F., Ahmed, S., Dey, S. S., & Wolsey, L. A., “Covering linear programming with violations,” INFORMS Journal on Computing, vol26, no.3, pp.531-546, 2014.
[2]	Porras, Á., Domínguez, C., Morales, J. M., & Pineda, S., “Tight and compact sample average approximation for joint chance constrained optimal power flow,” arXiv preprint arXiv:2205.03370, 2022.
[3]	Z. Wang, C. Shen, F. Liu, et al. Chance-constrained economic dispatch with non-Gaussian correlated wind power uncertainty. IEEE Trans. Power Syst., vol. 34, no. 3, pp. 4480-4893, 2017.
