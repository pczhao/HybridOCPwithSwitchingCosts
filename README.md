# HybridOCPwithSwitchingCosts
MATLAB code for solving hybrid optimal control problems with switching costs and free initial condition using occupation measures.

Requires SPOTless (found at <https://github.com/spot-toolbox/spotless>) and MOSEK (found at <https://www.mosek.com>).

To run the examples, run the run\_\*.m file.

To perform control synthesis, replace the spotsosprog.m file (located at spotless/spotopt/@spotsosprog) with the one provided (hybridOCP/src/tools/spotsosprog_changes.m). Or, make the following chages in the spotsosprog.m file by hand:
* Line 413: 
old: sol = minimize(pr,varargin)
new: [sol, y, basis, dual_multiplier] = minimize(pr,varargin)
* save all of the eqMultFac values in an array (look at lines 449, 454, 464, 471, 478 in the included spotsosprog_changes.m)
* Lines 84-92:  buildSOSDecompPrimal function
comment out lines 84-92


The code was tested with MATLAB R2018a and MOSEK 8.1.0.63 on an Intel Core i9, 18 core, 2.60 GHz, 128 GB RAM machine running Ubuntu 18.04 LTS.

* To solve hybrid ocp using occupation measure, run examples/run_T\*\_freeIC\_deg\*.m in MATLAB

* To solve hybrid ocp using GPOPS-II with **given** sequence of transition, run examples/Gpops/run\_Gpops\_T\*\_freeIC\_Main.m in MATLAB.
This code does not search through all possible sequences of transitions and therefore the computation time will be significantly smaller.

* To plot the trajectories, run examples/Plot/GenFigures.m in MATLAB.
