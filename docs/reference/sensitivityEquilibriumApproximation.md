# Equilibrium state approximation of the solution sensitivity for ODE systems

In this context, the sensitivity S(t;x,p) is dx(t;p)/dp, where x(t;p) is
the parameterized solution to an initial value problem for ordinary
differential equations and t is the independent varibale: x'=f(t,x;p),
where «'» indicates the derivative with respect to t. In cases where you
have a proxy variable for p, e.g. r=log(p), the chain rule applies.
Similarly, we also have an output sensitivity for the function
g(x(t;p)). The equilibrium approximation is exact for state-variable
values close to an equilibrium point q(p) (fixed-point): f(t,q(p);p)=0.

## Usage

``` r
sensitivityEquilibriumApproximation(
  experiments,
  model,
  parMap = identity,
  parMapJac = 1
)
```

## Arguments

- experiments:

  a list of simulation experiments

- model:

  a list of functions for the model the experiments are applicable to

- parMap:

  a map to transform parMCMC into p, parameters the model accepts

- simulations:

  an equivalent list of simulation results, for one parameter vector

- parMCMC:

  the parameters that are used in Markov chain Monte Carlo as the MC
  variable

## Value

a function S(parMCMC) -\> simulations_with_sensitivity, which attaches
the state sensitivity matrix array length(x) × length(p) × length(t) to
the simulations (solutions to the ODE).

## Details

The state sensitivity matrix:

               d state(time[k],state, param)[i]
    S[i,j,k] = --------------------------------  ,
               d param[j]

where param are the raw model parameters. This matrix is calculated as
an intermediate and then transformed into:

                d func(time[k], state, c(parMap(parMCMC),input))[i]
    Sh[i,j,k] = --------------------------------------------------
                d parMCMC[j]

where parMCMC is the Markov chain variable and usually shorter than
param as we typically don't sample all of the model's parameters. Some
model parameters may be known, some may be input parameters not
intrinsic to the model but related to the experimental setup (that is
why parMCMC and param are different).

This transformation requires the output function jacobian (funcJac) and
the parameter jacobian (funcJacp) in the model variable.

As we transform the parameters themselves, the chain rule requests
parMapJacl,k = d paraml / d parMCMCk

Typically, the sensitivity needs to be known at different time-points
t_k. The 3-dimensional array Si,j,k, where the index k corrsponds to
time t_k; the closer x(t_k) is to equilibrium, the better the
approximation; near the initial state, the sensitivity is also correct
(only the intermediate time-span is approximate).

This function requires pracma::expm to work.

The f

## Examples

``` r
if (FALSE) {
y <- simulate(parMCMC)
S <- sensitivityEquilibriumApproximation(experiments, model, parMap, parMapJac)
y <- S(parMap,y)
}
```
