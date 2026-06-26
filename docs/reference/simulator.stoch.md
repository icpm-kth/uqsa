# Function that creates a closure that simulates a stochastic trajectory with the Gillespie algorithm given certain experimental conditions and a parameter vector, and computes the distance between the simulation and the experimental data

Function that creates a closure that simulates a stochastic trajectory
with the Gillespie algorithm given certain experimental conditions and a
parameter vector, and computes the distance between the simulation and
the experimental data

## Usage

``` r
simulator.stoch(
  experiments,
  model.tab = model.tab,
  reactions = NULL,
  parMap = identity,
  outputFunction = function(t, state, param) {
     state
 },
  vol = 4e-16,
  unit = 1e-06,
  nStochSim = 3,
  distance = NULL
)
```

## Arguments

- experiments:

  a list of experiment.

- model.tab:

  an SBtab model.

- reactions:

  a list of GillespieSSA2 reactions.

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

- nStochSim:

  number of stochastic simulations to average over.

- distance:

  (optional) function that computes the distance between experimental
  data and simulated data. Default value is NULL.

## Value

a function that, given a parameter, returns a simulated trajectory
obtained via the Gillespie algorithm. Specifically, experimental time
points, output value measured at these time points, and (if distance is
not NULL) the distance between the simulated trajectories and the
experimental data.
