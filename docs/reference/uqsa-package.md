# uqsa: Uncertainty Quantification and Global Sensitivity Analysis

In the field of systems biology, chemical reaction networks are modeled
in various ways, two of those are: (i) stochastic simulations (e.g.
Gillespie algorithm) and (ii) ordinary differential equations. In this
package we use a simple tabular model description of reaction systems
and automatically generate C code for either solver type. We use the ODE
solvers from the GNU Scientific Library and provide an interface that
deals with lists of simulation experiments. Each simulation experiment
contains both the data, and instructions for the model to replicate the
data. We use ABC methods (combined with MCMC and SMC) as well as classic
MCMC methods such as Random Walk Metropolis (Gaussian transition kernel)
and SMMALA (Simplified Manifold Metropolis adjusted Langevin algorithm)
for a Bayesian investigation of the model´s parameter space. Experiments
can be evaluated in a sequence; intermediate probability densities are
modeled using the VineCopula package. UQSA is also intended to be useful
in an HPC environment, with some functions that use pbdMPI capabilities.

## See also

Useful links:

- <https://icpm-kth.github.io/uqsa/>

## Author

**Maintainer**: Andrei Kramer <andreikr@kth.se>

Authors:

- Alexandra Jauhiainen <Alexandra.Jauhiainen@astrazeneca.com>

- Olivia Eriksson <olivia@kth.se> \[contributor\]

- Federica Milinanni <fedmil@kth.se>
