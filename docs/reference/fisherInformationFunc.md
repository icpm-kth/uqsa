# Fisher Information from Sensitivity

Given a list of simulation sensitivities, this function returns the
fisher information (sum over all experiments). The actual work is done
in the returned function that implicitly depends on the model,
experiments, and parameter mapping

## Usage

``` r
fisherInformationFunc(
  experiments,
  parMap = identity,
  parMapJac = function(x) {
     diag(1, length(x))
 }
)
```

## Arguments

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
