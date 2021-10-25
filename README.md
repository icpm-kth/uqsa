# uqsa
Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

The code distributed in this repository implements the methodology presented in the paper "Uncertainty quantification, propagation and
characterization by Bayesian analysis combined with global sensitivity analysis applied to dynamical intracellular pathway models" by Eriksson & Jauhiainen et al (2018). The code is distributed under the GNU General Public License v3.0.

The UQ folder contains R scripts to run the uncertainty quantification method (ABC-MCMC with copulas). The packages ks, VineCopula, MASS, R.utils, and R.matlab are required (the last package to save the output data also into a MATLAB compatible file). The main script to run is called runABCMCMC-Phenotype123.R. This script will fit the model to phenotypes 1-3 (as described in the paper), which we use as an illustrative test case. The resulting data is also uploaded in the folder, both in R and MATLAB format. We use phenotype 4 as our prediction dataset to illustrate the SA methodology. 

Th SA folder contains MATLAB scripts to run global sensitivity analysis. A version of MATLAB later than 2014a is required. The main script to run is called get_predictions_do_SA.m and requires access to the file Draws-Phenotype123-Scale1000.mat which is available in the UQ folder. 

### Import an SBtab file in R

To import an SBtab (collection of) file(s) it is possible to use the function ´import_from_SBtab´ defined in the file "import_from_SBtab.R".

For information about the SBtab data format please refer to [the official git repository](https://github.com/tlubitz/SBtab).

Input to function ´ímport_from_SBtab´:
* the directory ´SBtabDir´ in which the SBtab spreadsheet are saved as tsv files.

Output of function ´import_from_SBtab´:
* an R variable with information from the SBtab spreadsheets

By running function ´ímport_from_SBtab´, a ´.vf´ is created in the same folder given as input to the function (´SBtabDir´).

Example
´sbtab_var <- import_from_SBtab(SBtabDir)´

### Full ODE-model

The full ODE-model correponding to all phenotypes can be recieved from the correponding authors upon request.

### Data references
The experimental data is extracted from the following references:<br/>
Stemmer PM, Klee CB. Biochemistry. 1994;33(22):6859-6866 (phenotype 1, 3 and 4)<br/>
O'Donnell SE et al. Proteins. 2011;79(3):765-786 (phenotype 2)

### Acknowledgements
This open source software code was developed in part or in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1 and SGA2).
