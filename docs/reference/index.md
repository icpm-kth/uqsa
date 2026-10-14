# Package index

## Model and data

Functions for loading model and data and construct a mathematical model
(ODE or stochastic) and corresponding code.

- [`uqsa_example()`](https://icpm-kth.github.io/uqsa/reference/uqsa_example.md)
  : Loads an example model for this package
- [`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)
  : Loads the content from a series of TSV files
- [`as_ode()`](https://icpm-kth.github.io/uqsa/reference/as_ode.md) :
  Interprets a model as an ODE
- [`as_cme()`](https://icpm-kth.github.io/uqsa/reference/as_cme.md) :
  Interprets the provided biological reaction model as stochastic
- [`write_and_compile()`](https://icpm-kth.github.io/uqsa/reference/write_and_compile.md)
  : Writes code to a file and compiles
- [`values()`](https://icpm-kth.github.io/uqsa/reference/values.md) :
  Finds the values in a data.frame
- [`uncertainty()`](https://icpm-kth.github.io/uqsa/reference/uncertainty.md)
  : Finds the uncertainty of values in a data.frame
- [`experiments()`](https://icpm-kth.github.io/uqsa/reference/experiments.md)
  : Extracts measured data and simulation experiment instructions
- [`standard_error_matrix()`](https://icpm-kth.github.io/uqsa/reference/standard_error_matrix.md)
  : Standard Error Matrix from an errors object
- [`conservation_law_analysis()`](https://icpm-kth.github.io/uqsa/reference/conservation_law_analysis.md)
  : Reduces the size of the system
- [`kinetic_law_matrix()`](https://icpm-kth.github.io/uqsa/reference/kinetic_law_matrix.md)
  : Represents the kinetic laws of the model as a two column matrix
- [`replace_powers()`](https://icpm-kth.github.io/uqsa/reference/replace_powers.md)
  : Interprets caret notation of arithmetic powers for C
- [`stoichiometry()`](https://icpm-kth.github.io/uqsa/reference/stoichiometry.md)
  : Returns a list of named stoichiometric vectors
- [`stoichiometric_matrix()`](https://icpm-kth.github.io/uqsa/reference/stoichiometric_matrix.md)
  : Infers the stoichiometric matrix of a reaction network
- [`unit.from.string()`](https://icpm-kth.github.io/uqsa/reference/unit.from.string.md)
  : Interprets a unit string
- [`unit.id()`](https://icpm-kth.github.io/uqsa/reference/unit.id.md) :
  Converts a unit to a string that works as an identifier
- [`unit_as_character()`](https://icpm-kth.github.io/uqsa/reference/unit_as_character.md)
  : Converts a unit data.frame into a printable string
- [`units_from_table()`](https://icpm-kth.github.io/uqsa/reference/units_from_table.md)
  : Finds units column in a data.frame and returns the values
- [`yJacobian()`](https://icpm-kth.github.io/uqsa/reference/yJacobian.md)
  : Calculates the Jacobian of string-math
- [`` `c_path<-`() ``](https://icpm-kth.github.io/uqsa/reference/c_path-set.md)
  : Adds information about the model's C code
- [`c_path()`](https://icpm-kth.github.io/uqsa/reference/c_path.md) :
  Retrieves information about the model's C code
- [`generate_code()`](https://icpm-kth.github.io/uqsa/reference/generate_code.md)
  : Constructs Code
- [`print(`*`<cme>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.cme.md)
  : Prints a Summary about the CME model
- [`print(`*`<ode>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.ode.md)
  : Prints a summary about the ode
- [`shlib()`](https://icpm-kth.github.io/uqsa/reference/shlib.md) :
  Compiles C code to shared library
- [`` `so_path<-`() ``](https://icpm-kth.github.io/uqsa/reference/so_path-set.md)
  : Adds information about compiled code
- [`so_path()`](https://icpm-kth.github.io/uqsa/reference/so_path.md) :
  Retrieves information about compiled code
- [`write_c_code()`](https://icpm-kth.github.io/uqsa/reference/write_c_code.md)
  : Writes C code to a file

## Simulations

Functions for simulating the model

- [`simulator.c()`](https://icpm-kth.github.io/uqsa/reference/simulator.c.md)
  : Creates a closure that simulates the model
- [`simstoch()`](https://icpm-kth.github.io/uqsa/reference/simstoch.md)
  : Simulates a Stochastic Model
- [`simfi()`](https://icpm-kth.github.io/uqsa/reference/simfi.md) :
  Creates a closure that simulates the model, similar to simulator.c
- [`scrnn()`](https://icpm-kth.github.io/uqsa/reference/scrnn.md) :
  Returns a closure around gsl_odeiv2_CRNN()
- [`gsl_odeiv2_fi()`](https://icpm-kth.github.io/uqsa/reference/gsl_odeiv2_fi.md)
  : Simulates an ode model
- [`gsl_odeiv2_CRNN()`](https://icpm-kth.github.io/uqsa/reference/gsl_odeiv2_CRNN.md)
  : Simulates a CRNN ode model
- [`name_method()`](https://icpm-kth.github.io/uqsa/reference/name_method.md)
  : Reverse look-up of method name from key

## Prior

Functions for constructing priors

- [`rNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/rNormalPrior.md)
  : Returns a random vector generator
- [`rUniformPrior()`](https://icpm-kth.github.io/uqsa/reference/rUniformPrior.md)
  : Returns a random vector generator
- [`rCopulaPrior()`](https://icpm-kth.github.io/uqsa/reference/rCopulaPrior.md)
  : Returns a function that generates random values from the copula
  model
- [`dNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/dNormalPrior.md)
  : Creates the density function of a multivariate normal distribution
  with independent components
- [`dUniformPrior()`](https://icpm-kth.github.io/uqsa/reference/dUniformPrior.md)
  : Creates a uniform density function
- [`dCopulaPrior()`](https://icpm-kth.github.io/uqsa/reference/dCopulaPrior.md)
  : Creates a prior probability density function
- [`gNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/gNormalPrior.md)
  : Creates the gradient function of a multivariate normal distribution
  with independent components, in log-space
- [`fitCopula()`](https://icpm-kth.github.io/uqsa/reference/fitCopula.md)
  : Makes a Probability Density Estimator from a Sample

## Uncertainty quantification - likelihood based

MCMC methods for Bayesian parameter inference

- [`mcmc()`](https://icpm-kth.github.io/uqsa/reference/mcmc.md) :
  Creates a Markov Chain Monte Carlo sampler
- [`mcmc_init()`](https://icpm-kth.github.io/uqsa/reference/mcmc_init.md)
  : Initializes the Markov chain
- [`mcmc_mpi()`](https://icpm-kth.github.io/uqsa/reference/mcmc_mpi.md)
  : Creates the MPI version of an MCMC sampler
- [`metropolis_update()`](https://icpm-kth.github.io/uqsa/reference/metropolis_update.md)
  : Performs an MCMC update using the Metropolis rule
- [`high_level_metropolis()`](https://icpm-kth.github.io/uqsa/reference/high_level_metropolis.md)
  : Determines a likely Metropolis setup based on the information in the
  model
- [`smmala_update()`](https://icpm-kth.github.io/uqsa/reference/smmala_update.md)
  : Performs an SMMALA update
- [`high_level_smmala()`](https://icpm-kth.github.io/uqsa/reference/high_level_smmala.md)
  : Determines a likely SMMALA setup based on the information in the
  model
- [`logLikelihoodFunc()`](https://icpm-kth.github.io/uqsa/reference/logLikelihoodFunc.md)
  : Constructs a custom log-likelihood function
- [`ll()`](https://icpm-kth.github.io/uqsa/reference/ll.md) : Default
  log-likelihood function
- [`gllf()`](https://icpm-kth.github.io/uqsa/reference/gllf.md) :
  Default gradient-log-likelihood function
- [`fi()`](https://icpm-kth.github.io/uqsa/reference/fi.md) : Default
  Fisher-information function
- [`tune_step_size()`](https://icpm-kth.github.io/uqsa/reference/tune_step_size.md)
  : Finds a good step-size for a given MCMC algorithm
- [`loadSample_mpi()`](https://icpm-kth.github.io/uqsa/reference/loadSample_mpi.md)
  : Merges mpi-samples into one
- [`gatherSample()`](https://icpm-kth.github.io/uqsa/reference/gatherSample.md)
  : Collects all sample points, from all files, with the given
  temperature
- [`gatherReplicas()`](https://icpm-kth.github.io/uqsa/reference/gatherReplicas.md)
  : Collects statistical Replicas

## Uncertainty quantification - ABC

Approximate Bayesian Computation (ABC) MCMC for parameter inference

- [`makeObjective()`](https://icpm-kth.github.io/uqsa/reference/makeObjective.md)
  : Creates objective functions from ingredients
- [`abc_mcmc()`](https://icpm-kth.github.io/uqsa/reference/abc_mcmc.md)
  : Performs Approximate Bayesian Computation Sampling of Model
  Parameters via Markov Chain Monte Carlo
- [`ABCSMC()`](https://icpm-kth.github.io/uqsa/reference/ABCSMC.md) :
  Performs Approximate Bayesian Computation as a Particle Filter

## Global Sensitivity Analysis

- [`gsa_binning()`](https://icpm-kth.github.io/uqsa/reference/gsa_binning.md)
  : Performs global sensitivity analysis using the binning method
- [`gsa_saltelli()`](https://icpm-kth.github.io/uqsa/reference/gsa_saltelli.md)
  : Outputs the global sensitivity scores SI and SIT, calculated by the
  Sobol-Homma-Saltelli method
- [`saltelli_prior()`](https://icpm-kth.github.io/uqsa/reference/saltelli_prior.md)
  : Creates prior samples for the Sobol-Homma-Saltelli global
  sensitivity analysis
- [`KLD()`](https://icpm-kth.github.io/uqsa/reference/KLD.md) : Kullback
  Leibler Divergence
- [`exact_normal_kld()`](https://icpm-kth.github.io/uqsa/reference/exact_normal_kld.md)
  : Multivariate Normal Distribution KLD

## Parameter mappings

Parameter transformations from sampling space back to model space

- [`logParMap()`](https://icpm-kth.github.io/uqsa/reference/logParMap.md)
  : NATURAL LOG parameter mapping used by the MCMC module
- [`logParMapJac()`](https://icpm-kth.github.io/uqsa/reference/logParMapJac.md)
  : NATURAL LOG parameter mapping, jacobian
- [`log10ParMap()`](https://icpm-kth.github.io/uqsa/reference/log10ParMap.md)
  : LOG10 parameter mapping used by the MCMC module
- [`log10ParMapJac()`](https://icpm-kth.github.io/uqsa/reference/log10ParMapJac.md)
  : LOG10 parameter mapping, Jacobian
- [`log2ParMap()`](https://icpm-kth.github.io/uqsa/reference/log2ParMap.md)
  : LOG2 parameter mapping used by the MCMC module
- [`log2ParMapJac()`](https://icpm-kth.github.io/uqsa/reference/log2ParMapJac.md)
  : LOG2 parameter mapping, Jacobian

## Helper Functions and Operators

- [`` `%@%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-at-grapes.md)
  : Fetches an Attribute, with checks
- [`` `%as%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-as-grapes.md)
  : A binary operator on strings with units in them
- [`` `%has%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-has-grapes.md)
  : Checks whether a variable has the named attributes
- [`` `%otherwise%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-otherwise-grapes.md)
  : Specifies default values, as fallback
- [`determinePrefix()`](https://icpm-kth.github.io/uqsa/reference/determinePrefix.md)
  : Determines a prefix from a character vector str of similar contents
- [`` `modify<-`() ``](https://icpm-kth.github.io/uqsa/reference/modify-set.md)
  : Modifies a value
- [`` `base<-`() ``](https://icpm-kth.github.io/uqsa/reference/base-set.md)
  : Converts values to linear space
- [`column()`](https://icpm-kth.github.io/uqsa/reference/column.md) :
  Retrieves j-th column with names
- [`print(`*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.experiments.md)
  : Prints the simulation experiments
- [`print(`*`<mcmcVariable>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.mcmcVariable.md)
  : Prints information about the mcmc variable
- [`print(`*`<simulation>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.simulation.md)
  : Prints the simulation results
- [`print(`*`<unit_of_measurement>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.unit_of_measurement.md)
  : Prints an interpretation string of a unit
- [`small()`](https://icpm-kth.github.io/uqsa/reference/small.md) :
  Reduces the sample to its effective size

## Plotting functions

- [`highCor()`](https://icpm-kth.github.io/uqsa/reference/highCor.md) :
  Orders index-pairs from high to low correlation
- [`pc_dist()`](https://icpm-kth.github.io/uqsa/reference/pc_dist.md) :
  Plots a sample in parallel coordinates
- [`sensitivity.graph()`](https://icpm-kth.github.io/uqsa/reference/sensitivity.graph.md)
  : Plot the result of global sensitivity analysis
- [`showPosterior()`](https://icpm-kth.github.io/uqsa/reference/showPosterior.md)
  : Makes a pairs plot for a sample using density estimation methods

## Internal functions

- [`onlyCoefficients()`](https://icpm-kth.github.io/uqsa/reference/onlyCoefficients.md)
  : Returns a list of reaction coefficients
- [`onlyNames()`](https://icpm-kth.github.io/uqsa/reference/onlyNames.md)
  : Returns only the names in a reaction formula
- [`CRNN()`](https://icpm-kth.github.io/uqsa/reference/CRNN.md) :
  Creates C code for a chemical reaction neural network
- [`clear_yacas_environment()`](https://icpm-kth.github.io/uqsa/reference/clear_yacas_environment.md)
  : Clears Yacas variables
- [`defaultDistance()`](https://icpm-kth.github.io/uqsa/reference/defaultDistance.md)
  : Default distance function for one experiment
- [`formulae()`](https://icpm-kth.github.io/uqsa/reference/formulae.md)
  : Finds a column that contains some kind of mathematic expression in a
  data.frame
- [`linear_scale()`](https://icpm-kth.github.io/uqsa/reference/linear_scale.md)
  : Interprets a character vector as names of logarithms
- [`parse_concise()`](https://icpm-kth.github.io/uqsa/reference/parse_concise.md)
  : Reads Concise Error Notation
- [`method()`](https://icpm-kth.github.io/uqsa/reference/method.md) :
  Finds the method id given a named method
- [`plot(`*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/plot.experiments.md)
  : Plot function for experiments
- [`` `[`( ``*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/sub-.experiments.md)
  : Subsets experiments with preserved class
- [`` `[`( ``*`<simulation>`*`)`](https://icpm-kth.github.io/uqsa/reference/sub-.simulation.md)
  : Subset simulations with preserved class
