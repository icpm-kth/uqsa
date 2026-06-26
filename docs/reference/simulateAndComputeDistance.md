# Function that simulates a stochastic trajectory with the Gillespie algorithm given certain experimental conditions and a parameter vector, and computes the distance between the simulation and the experimental data

Function that simulates a stochastic trajectory with the Gillespie
algorithm given certain experimental conditions and a parameter vector,
and computes the distance between the simulation and the experimental
data

## Usage

``` r
simulateAndComputeDistance(
  e,
  param,
  parMap = parMap,
  Phi = Phi,
  parameters_from_expressions = parameters_from_expressions,
  nStochSim = nStochSim,
  reactions = compiled_reactions,
  distance = distance
)
```

## Arguments

- e:

  an experiment

- param:

  a named parameter vector

- parMap:

  a function that translates ABC variables (parABC) into something the
  model will accept.

- Phi:

  Volume

- parameters_from_expressions:

  a vector of evaluated expressions

- nStochSim:

  number of stochastic simulations to average over

- reactions:

  a list that encodes the reactions for GillespieSSA2

- distance:

  a user supplied function that calculates a distance between simulation
  and data with an interface of distance(simulation, data, errVal),
  where errVal is an estimate of the measuremnet noise (e.g. standard
  deviation), if needed by the function.

## Value

the distance between the trajectory just simulated and the experiment
considered
