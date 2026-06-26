# Function that creates the objective function

Given a parameter set, this function computes the distance between
experimental data and simulated data (coresponding to the parameter in
input).

## Usage

``` r
makeObjectiveSSA(
  experiments,
  model.tab,
  parNames,
  distance,
  parMap = identity,
  outputFunction = identity,
  vol = 4e-16,
  unit = 1e-06,
  reactions,
  nStochSim = 1
)
```

## Arguments

- experiments:

  a list of experiments.

- model.tab:

  an SBtab model.

- parNames:

  the names of the (biological) parameters of the model.

- distance:

  a user supplied function that calculates a distance between simulation
  and data with an interface of distance(simulation, data, errVal),
  where errVal is an estimate of the measurement noise (e.g. standard
  deviation), if needed by the function.

- parMap:

  a function that translates ABC variables (parABC) into something the
  model will accept.

- outputFunction:

  a function that, given (t,state,param), outputs the quantity that is
  measured in the experiments (e.g., the outputFunction may output the
  concentration of one of the compounds in the system).

- vol:

  Volume in which the reactions take place.

- unit:

  unit of measure for the volume (it converts the volume into liters).

- reactions:

  a list of GillespieSSA2 reactions

- nStochSim:

  number of stochastic simulations to average over.

## Value

a closure for the objective function that implicitly depends on all of
the arguments to this function but explicitly only on the ABC parameters
parABC.
