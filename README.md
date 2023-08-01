# Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

This is an R package that performs *parameter estimation*,
*uncertainty quantification*, and *global sensitivity analysis* using
Bayesian methods and ABC-MCMC sampling.

## Installation

```R
# requires the 'remotes' package
remotes::install_github("icpm-kth/uqsa")
```

Alternatively, download a release `tar.gz` or `zip` file and run

```sh
R CMD INSTALL uqsa*.{tar.gz,zip}
```

## Documentation

A detailed documentation is available in the [wiki](https://github.com/icpm-kth/uqsa/wiki/Documentation).

## Demos

A simple example using the AKAR4 model is available as a Jupyter notebook in the file [`UQ_AKAR4.ipynb`](https://github.com/icpm-kth/uqsa/blob/master/UQ_AKAR4.ipynb). You can run it locally after installing Jupyter (read [here](https://jupyter.org/install) for further imformation).

Three larger examples are available in the form of R scripts and can be run locally or (preferably) on a computer cluster:
* [AKAP79 ODE model](https://github.com/icpm-kth/uqsa/blob/master/inst/extdata/AKAP79/runABCMCMC_AKAP79.R)
* [CaMKII ODE model](https://github.com/icpm-kth/uqsa/blob/master/inst/extdata/CaMKII/runABCMCMC_CaMKII.R)
* [AKAR4 stochastic model](https://github.com/icpm-kth/uqsa/blob/master/inst/extdata/AKAR4/runABCMCMC_AKAR4_withStochasticStimulation.R)
  
These R scrips can be found in the folder [`inst/extdata`](https://github.com/icpm-kth/uqsa/tree/master/inst/extdata).
  
## Original Work

This repository (https://github.com/icpm-kth/uqsa) is a continued development of the methodology presented in the paper:
```
Eriksson, Olivia and Jauhiainen Alexandra et al. "Uncertainty quantification, propagation and characterization by 
Bayesian analysis combined with global sensitivity analysis applied to dynamical intracellular pathway models." 
Bioinformatics 35.2 (2019): 284-292.
```
The original code for that paper can be found at: https://github.com/alexjau/uqsa.

The demo UQ.ipynb was presented in Figure 5c of the paper
```
Eriksson, Olivia et al. "Combining hypothesis- and data-driven neuroscience modeling in FAIR workflows", 
eLife (2022), 11:e69013 DOI: https://doi.org/10.7554/eLife.69013
```
 The code is distributed under the GNU General Public License v3.0.

## Data references

The experimental data for the AKAP79 and AKAR4 models are from the publication:
```
Church, Timothy W., et al. "AKAP79 enables calcineurin to directly suppress protein kinase A activity." 
Elife 10 (2021): e68164.
```
and described in detail in the SBtab files of the AKAP79 and AKAR4 model folders.

The experimental data for the CaMKII model are extracted from the following references:
```
Stemmer PM, Klee CB. Biochemistry. 1994;33(22):6859-6866 (phenotype 1, 3 and 4)
```
```
O'Donnell SE et al. Proteins. 2011;79(3):765-786 (phenotype 2)
```


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

## Acknowledgements

This open source software code was developed in part or in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1 and SGA2).
