# Fisher Information from Sensitivity

Given a list of simulation sensitivities, this function returns the
fisher information (sum over all experiments). The actual work is done
in the returned function that implicitly depends on the model,
experiments, and parameter mapping

## Usage

``` r
fisherInformation(model, experiments, parMap = identity, parMapJac = 1)
```

## Arguments

- model:

  list of R functions for the ODE model

- experiments:

  list of experiments, with inputs

- parMap:

  mapping between MCMC variables and ODE parameters

- parMapJac:

  the jacobian of the above map

## Value

fisher information calculating funciton

## Details

return value: function(par, simulations, sensitivity) -\>
fisherInformation (matrix)

where par refers to the model parameters (possibly transformed), and
simulations performed with those parameters.
