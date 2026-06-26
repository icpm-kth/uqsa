# Determine a starting value for ABC's delta

In ABC settings a model solution is compared to data with an acceptance
threshold: delta. This-pre calibration function attempts to adjust this
delta value.

## Usage

``` r
preCalibration(
  objectiveFunction,
  npc = 1000,
  rprior,
  rep = 1,
  p = 0.05,
  sfactor = 0.1,
  delta = 0.01,
  num = 1
)
```

## Arguments

- objectiveFunction:

  function that, given a (vectorial) parameter as input, (1) simulates
  the model with the given parameter, and (2) outputs the distance
  between experimental data and simulated data simulated.

- npc:

  sample size of pre-calibration.

- rprior:

  a function that generates random ABC variables, distributed according
  to the prior.

- rep:

  number of repetitions of the preCalibration process.

- p:

  fraction (top scoring) of sampled points to base Sigma on (Sigma is
  the covariance matrix for the moves proposed in the ABCMCMC
  algorithm).

- sfactor:

  scales Sigma up or down (Sigma is the covariance matrix for the moves
  proposed in the ABCMCMC algorithm).

- delta:

  ABC threshold.

- num:

  number of different starting parameter vectors (initial states of the
  chains) to generate. Usually, num is equal to the number of chain that
  will be run in the sampling procedure.

## Value

list with entries prePar (sampled parameters), preDelta (distances
between experimental data and trajectories produced with each of the
parameters in prePar), Sigma (covariance matrix for the moves proposed
in the ABCMCMC algortihm) and startPar (starting parameters for the
ABCMCMC chains)

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
cme <- as_cme(m)
ex <- experiments(m)
C <- generate_code(cme)
c_path(cme) <- write_c_code(C)
so_path(cme) <- shlib(cme)
p0 <- log10(values(m$Parameter))
s <- simstoch(ex,cme,parMap=log10ParMap)
Obj <- makeObjective(ex,s)
rprior <- rNormalPrior(p0,0.5)
PC <- preCalibration(Obj,50,rprior,delta=1) # should be more than 50
#> Warning: distances between experiment and simulation are too big; selecting the best (50) parameter vectors.
print(names(PC))
#> [1] "prePar"   "preDelta" "Sigma"    "startPar"
print(PC$Sigma)
#>             [,1]        [,2]        [,3]
#> [1,] 0.004496082 0.001051894 0.003292040
#> [2,] 0.001051894 0.024296263 0.005356373
#> [3,] 0.003292040 0.005356373 0.029494625
```
