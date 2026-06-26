# Package index

## Model and data

Functions for loading model and data and construct a mathematical model
(ODE or stochastic) and corresponding code.

- [`uqsa_example()`](https://icpm-kth.github.io/uqsa/reference/uqsa_example.md)
  : Load an example model for this package
- [`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)
  : model_from_tsv loads the content from a series of tsv files
- [`as_ode()`](https://icpm-kth.github.io/uqsa/reference/as_ode.md) :
  Interpret a model as an ODE
- [`as_cme()`](https://icpm-kth.github.io/uqsa/reference/as_cme.md) :
  Interprets the provided model as a stochastic model
- [`write_and_compile()`](https://icpm-kth.github.io/uqsa/reference/write_and_compile.md)
  : Writes code to file and compiles
- [`values()`](https://icpm-kth.github.io/uqsa/reference/values.md) :
  Find values in a data.frame that is derived from a tsv file or similar
- [`uncertainty()`](https://icpm-kth.github.io/uqsa/reference/uncertainty.md)
  : Find the uncertainty of values in a data.frame that is derived from
  a tsv file or similar
- [`experiments()`](https://icpm-kth.github.io/uqsa/reference/experiments.md)
  : Extract Measured Data and Simulation Experiment Instructions
- [`standard_error_matrix()`](https://icpm-kth.github.io/uqsa/reference/standard_error_matrix.md)
  : Standard Error Matrix from an errors object
- [`conservation_law_analysis()`](https://icpm-kth.github.io/uqsa/reference/conservation_law_analysis.md)
  : Reduce the size of the system
- [`kinetic_law_matrix()`](https://icpm-kth.github.io/uqsa/reference/kinetic_law_matrix.md)
  : Split Kinetic Law
- [`replace_powers()`](https://icpm-kth.github.io/uqsa/reference/replace_powers.md)
  : replace_powers
- [`stoichiometry()`](https://icpm-kth.github.io/uqsa/reference/stoichiometry.md)
  : This function returns a list of named stoichiometric vectors
- [`stoichiometric_matrix()`](https://icpm-kth.github.io/uqsa/reference/stoichiometric_matrix.md)
  : The stoichiometric matrix of a reaction network
- [`unit.from.string()`](https://icpm-kth.github.io/uqsa/reference/unit.from.string.md)
  : Unit Interpreter
- [`unit.id()`](https://icpm-kth.github.io/uqsa/reference/unit.id.md) :
  Converts a unit to a string that works as an identifier
- [`unit.info()`](https://icpm-kth.github.io/uqsa/reference/unit.info.md)
  : Prints an interpretation string of a unit
- [`unit_as_character()`](https://icpm-kth.github.io/uqsa/reference/unit_as_character.md)
  : converts a unit data.frame into a printable string
- [`units_from_table()`](https://icpm-kth.github.io/uqsa/reference/units_from_table.md)
  : Get units from a data.frame column
- [`yJacobian()`](https://icpm-kth.github.io/uqsa/reference/yJacobian.md)
  : Jacobian of string-math
- [`` `c_path<-`() ``](https://icpm-kth.github.io/uqsa/reference/c_path-set.md)
  : Add information about the model's C code
- [`c_path()`](https://icpm-kth.github.io/uqsa/reference/c_path.md) :
  Retrieve information about the model's C code
- [`generate_code()`](https://icpm-kth.github.io/uqsa/reference/generate_code.md)
  : Construct Code
- [`print(`*`<cme>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.cme.md)
  : Print a Summary about the CME model
- [`print(`*`<ode>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.ode.md)
  : Print a summary about the ode
- [`shlib()`](https://icpm-kth.github.io/uqsa/reference/shlib.md) :
  Compile C code to shared library
- [`` `so_path<-`() ``](https://icpm-kth.github.io/uqsa/reference/so_path-set.md)
  : Add information about compiled code
- [`so_path()`](https://icpm-kth.github.io/uqsa/reference/so_path.md) :
  Retrieve information about compiled code
- [`write_c_code()`](https://icpm-kth.github.io/uqsa/reference/write_c_code.md)
  : Write the C code to a file

## Simulations

Functions for simulating the model

- [`simulator.c()`](https://icpm-kth.github.io/uqsa/reference/simulator.c.md)
  : This creates a closure that simulates the model
- [`simstoch()`](https://icpm-kth.github.io/uqsa/reference/simstoch.md)
  : Simulate stochastic model
- [`simfi()`](https://icpm-kth.github.io/uqsa/reference/simfi.md) : This
  creates a closure that simulates the model, similar to simulator.c
- [`scrnn()`](https://icpm-kth.github.io/uqsa/reference/scrnn.md) :
  scrnn returns a closure around gsl_odeiv2_CRNN()
- [`gsl_odeiv2_fi()`](https://icpm-kth.github.io/uqsa/reference/gsl_odeiv2_fi.md)
  : simulates an ode model with extra work
- [`gsl_odeiv2_CRNN()`](https://icpm-kth.github.io/uqsa/reference/gsl_odeiv2_CRNN.md)
  : simulates a CRNN ode model with extra work
- [`name_method()`](https://icpm-kth.github.io/uqsa/reference/name_method.md)
  : Reverse look-up of method name from key

## Prior

Functions for constructing priors

- [`rNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/rNormalPrior.md)
  : rNormalPrior returns a random vector generator
- [`rUniformPrior()`](https://icpm-kth.github.io/uqsa/reference/rUniformPrior.md)
  : rUniformPrior returns a random vector generator
- [`rCopulaPrior()`](https://icpm-kth.github.io/uqsa/reference/rCopulaPrior.md)
  : rCopulaPrior returns a function that generates random values from
  the copula model
- [`dNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/dNormalPrior.md)
  : dNormalPrior creates the density function of a multivariate normal
  distribution with independent components
- [`dUniformPrior()`](https://icpm-kth.github.io/uqsa/reference/dUniformPrior.md)
  : dUniformPrior creates a uniform density function
- [`dCopulaPrior()`](https://icpm-kth.github.io/uqsa/reference/dCopulaPrior.md)
  : dCopulaPrior creates a prior probability density function
- [`gNormalPrior()`](https://icpm-kth.github.io/uqsa/reference/gNormalPrior.md)
  : gNormalPrior creates the gradient function of a multivariate normal
  distribution with independent components, in log-space
- [`fitCopula()`](https://icpm-kth.github.io/uqsa/reference/fitCopula.md)
  : Makes a Probability Density Estimate (from a sample)

## Uncertainty quantification - likelihood based

MCMC methods for Bayesian parameter inference

- [`mcmc()`](https://icpm-kth.github.io/uqsa/reference/mcmc.md) : Markov
  Chain Monte Carlo
- [`mcmc_init()`](https://icpm-kth.github.io/uqsa/reference/mcmc_init.md)
  : Initialize the Markov chain
- [`mcmc_mpi()`](https://icpm-kth.github.io/uqsa/reference/mcmc_mpi.md)
  : The MPI version of the mcmc function
- [`metropolis_update()`](https://icpm-kth.github.io/uqsa/reference/metropolis_update.md)
  : Metropolis Update is an MCMC update function
- [`high_level_metropolis()`](https://icpm-kth.github.io/uqsa/reference/high_level_metropolis.md)
  : High Level Metropolis function
- [`smmala_update()`](https://icpm-kth.github.io/uqsa/reference/smmala_update.md)
  : SMMALA Update is an MCMC update function
- [`high_level_smmala()`](https://icpm-kth.github.io/uqsa/reference/high_level_smmala.md)
  : High Level SMMALA function
- [`logLikelihoodFunc()`](https://icpm-kth.github.io/uqsa/reference/logLikelihoodFunc.md)
  : Default log-likelihood function
- [`ll()`](https://icpm-kth.github.io/uqsa/reference/ll.md) : Default
  Log-likelihood Function
- [`gllf()`](https://icpm-kth.github.io/uqsa/reference/gllf.md) :
  Default gradient-log-likelihood Function
- [`fi()`](https://icpm-kth.github.io/uqsa/reference/fi.md) : Default
  gradient-Log-likelihood Function
- [`tune_step_size()`](https://icpm-kth.github.io/uqsa/reference/tune_step_size.md)
  : Find a good Step-Size for a given MCMC Algorithm
- [`change_temperature()`](https://icpm-kth.github.io/uqsa/reference/change_temperature.md)
  : Should 2 Markov chains exchange their temperatures
- [`loadSample_mpi()`](https://icpm-kth.github.io/uqsa/reference/loadSample_mpi.md)
  : This function merges mpi-samples into one
- [`gatherSample()`](https://icpm-kth.github.io/uqsa/reference/gatherSample.md)
  : gatherSample collects all sample points, from all files, with the
  given temperature
- [`gatherReplicas()`](https://icpm-kth.github.io/uqsa/reference/gatherReplicas.md)
  : Collect statistical Replicas

## Uncertainty quantification - ABC

Approximate Bayesian Computation (ABC) MCMC for parameter inference

- [`makeObjective()`](https://icpm-kth.github.io/uqsa/reference/makeObjective.md)
  : creates Objective functions from ingredients
- [`preCalibration()`](https://icpm-kth.github.io/uqsa/reference/preCalibration.md)
  : Determine a starting value for ABC's delta
- [`ABCMCMC()`](https://icpm-kth.github.io/uqsa/reference/ABCMCMC.md) :
  Performs and Approximate Bayesian Computation Sampling of Model
  Parameters
- [`ABCSMC()`](https://icpm-kth.github.io/uqsa/reference/ABCSMC.md) :
  Performs and Approximate Bayesian Computation as a Particle Filter
- [`checkFitWithPreviousExperiments()`](https://icpm-kth.github.io/uqsa/reference/checkFitWithPreviousExperiments.md)
  : ABC acceptance of currently sampled values given old data (Prior)

## Global Sensitivity Analysis

- [`gsa_binning()`](https://icpm-kth.github.io/uqsa/reference/gsa_binning.md)
  : Global Sensitivity Analysis
- [`gsa_saltelli()`](https://icpm-kth.github.io/uqsa/reference/gsa_saltelli.md)
  : Outputs the global sensitivity scores SI and SIT, calculated by the
  Sobol-Homma-Saltelli method
- [`saltelli_prior()`](https://icpm-kth.github.io/uqsa/reference/saltelli_prior.md)
  : Sample for the Sobol-Homma-Saltelli Global Sensitivity Analysis

## Parameter mappings

Parameter transformations from sampling space back to model space

- [`logParMap()`](https://icpm-kth.github.io/uqsa/reference/logParMap.md)
  : NATURAL LOG parameter mapping used by the MCMC module
- [`logParMapJac()`](https://icpm-kth.github.io/uqsa/reference/logParMapJac.md)
  : NATURAL LOG parameter mapping, jacobian
- [`log10ParMap()`](https://icpm-kth.github.io/uqsa/reference/log10ParMap.md)
  : LOG10 parameter mapping used by the MCMC module
- [`log10ParMapJac()`](https://icpm-kth.github.io/uqsa/reference/log10ParMapJac.md)
  : LOG10 parameter mapping, jacobian
- [`log2ParMap()`](https://icpm-kth.github.io/uqsa/reference/log2ParMap.md)
  : LOG2 parameter mapping used by the MCMC module
- [`log2ParMapJac()`](https://icpm-kth.github.io/uqsa/reference/log2ParMapJac.md)
  : LOG2 parameter mapping, jacobian

## Helper Functions and Operators

- [`` `%@%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-at-grapes.md)
  : Fetch an Attribute
- [`` `%as%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-as-grapes.md)
  : %as% is a binary operator on strings with units in them
- [`` `%has%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-has-grapes.md)
  : checks whether a variable has the named attributes
- [`` `%otherwise%` ``](https://icpm-kth.github.io/uqsa/reference/grapes-otherwise-grapes.md)
  : This function can be used to specify default values
- [`determinePrefix()`](https://icpm-kth.github.io/uqsa/reference/determinePrefix.md)
  : Determine a prefix from a character vector str of similar contents
- [`` `modify<-`() ``](https://icpm-kth.github.io/uqsa/reference/modify-set.md)
  : Modifies a value
- [`` `base<-`() ``](https://icpm-kth.github.io/uqsa/reference/base-set.md)
  : Convert to linear space
- [`column()`](https://icpm-kth.github.io/uqsa/reference/column.md) :
  Get j-th column with names
- [`print(`*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.experiments.md)
  : prints the simulation experiments
- [`print(`*`<mcmcVariable>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.mcmcVariable.md)
  : print information about the mcmc variable
- [`print(`*`<simulation>`*`)`](https://icpm-kth.github.io/uqsa/reference/print.simulation.md)
  : prints the simulation results

## Plotting functions

- [`highCor()`](https://icpm-kth.github.io/uqsa/reference/highCor.md) :
  highCor returns ordered index-pairs of high to low correlation
- [`pcDist()`](https://icpm-kth.github.io/uqsa/reference/pcDist.md) :
  plots a sample in parallel coordinates
- [`sensitivity.graph()`](https://icpm-kth.github.io/uqsa/reference/sensitivity.graph.md)
  : plot the sensitivity matrix
- [`showPosterior()`](https://icpm-kth.github.io/uqsa/reference/showPosterior.md)
  : showPosterior makes a pairs plot for a sample

## Internal functions

- [`onlyCoefficients()`](https://icpm-kth.github.io/uqsa/reference/onlyCoefficients.md)
  : Returns a list of reaction coefficients

- [`onlyNames()`](https://icpm-kth.github.io/uqsa/reference/onlyNames.md)
  : Returns only the names in a reaction formula

- [`CRNN()`](https://icpm-kth.github.io/uqsa/reference/CRNN.md) : CRNN
  creates C code for a chemical reaction neural network

- [`clear_yacas_environment()`](https://icpm-kth.github.io/uqsa/reference/clear_yacas_environment.md)
  : Clear Yacas variables

- [`defaultDistance()`](https://icpm-kth.github.io/uqsa/reference/defaultDistance.md)
  : default distance function for one experiment

- [`formulae()`](https://icpm-kth.github.io/uqsa/reference/formulae.md)
  : Find a column that contains some kind of mathematic expression in a
  data.frame

- [`linear_scale()`](https://icpm-kth.github.io/uqsa/reference/linear_scale.md)
  : Interprets a character vector as names of logarithms

- [`parse_concise()`](https://icpm-kth.github.io/uqsa/reference/parse_concise.md)
  : Read Concise Error Notation

- [`` `reaction<-`() ``](https://icpm-kth.github.io/uqsa/reference/reaction-set.md)
  : Add a Reaction to an ODE

- [`update_values()`](https://icpm-kth.github.io/uqsa/reference/update_values.md)
  :

  Updates the named values of vector `v` with values mentioned in
  data.frame d

- [`simple.unit()`](https://icpm-kth.github.io/uqsa/reference/simple.unit.md)
  : Simple unit from string

- [`method()`](https://icpm-kth.github.io/uqsa/reference/method.md) :
  Find Integer

- [`plot(`*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/plot.experiments.md)
  : plot function for experiments

- [`` `[`( ``*`<experiments>`*`)`](https://icpm-kth.github.io/uqsa/reference/sub-.experiments.md)
  : Subset experiments with preserved class

- [`` `[`( ``*`<simulation>`*`)`](https://icpm-kth.github.io/uqsa/reference/sub-.simulation.md)
  : Subset simulations with preserved class
