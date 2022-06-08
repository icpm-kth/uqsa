# uqsa
Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

This repository (https://github.com/icpm-kth/uqsa) is a continued development of the methodology presented in the paper:
```
Eriksson, Olivia and Jauhiainen Alexandra et al. "Uncertainty quantification, propagation and characterization by 
Bayesian analysis combined with global sensitivity analysis applied to dynamical intracellular pathway models." 
Bioinformatics 35.2 (2019): 284-292.
```
The original code for that paper can be found at: https://github.com/alexjau/uqsa

The demo UQ.ipynb was presented in Figure 5c of the paper
```
Eriksson, Olivia et al. "Combining hypothesis- and data-driven neuroscience modeling in FAIR workflows", 
eLife (2022), 11:e69013 DOI: https://doi.org/10.7554/eLife.69013
```
 The code is distributed under the GNU General Public License v3.0.

The UQ folder contains R scripts to run the uncertainty quantification method (ABC-MCMC with copulas). A detailed documentation of this code can be found in the wiki of this repository (https://github.com/icpm-kth/uqsa/wiki/Documentation).

The SA folder will soon contain the code for the global sensitivity analysis, meanwhile we refer to https://github.com/alexjau/uqsa.

### Data references
The experimental data for the AKAP79 and AKAR4 models are from the publication:
```
Church, Timothy W., et al. "AKAP79 enables calcineurin to directly suppress protein kinase A activity." 
Elife 10 (2021): e68164.
```
and described in detail in the SBtab files of the AKAP79 and AKAR4 model folders.

### Model references
The AKAP79 model is from
```
Church, Timothy W., et al. "AKAP79 enables calcineurin to directly suppress protein kinase A activity." 
Elife 10 (2021): e68164.`

```
which was modified from:
```
Buxbaum JD, Dudai Y. 1989. A quantitative model for the kinetics of cAMP-dependent protein kinase (type II) 
activity. Long-term activation of the kinase and its possible relevance to learning and memory. 
The Journal of Biological Chemistry 264:9344–9351.
```

### Acknowledgements
This open source software code was developed in part within the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270, No. 785907 and No 945539 (Human Brain Project SGA1, SGA2 and SGA3); the Swedish Research Council VR-M-2017-02806 and VR-M-2020-01652; and Governmental grants for strategic research areas (Swedish e-Science Research Center; Digital Futures).
