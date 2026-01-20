# Uncertainty Quantification: run model ordinary differential equation
# Copyright (C) 2022 Federica Milinanni (fedmil@kth.se)

#' simulates an ode model with extra work
#'
#' This function calls a C function which solves an initial value
#' problem, calculates the sensitivity of the solution, log-likelihood
#' value `ll`, gradient of `ll` amd Fisher-Information.
#'
#' @param name the name of the ODE model to simulate (a shared library of the same name will be dynamically loaded and needs to be created first)
#' @param experiments a list of `N` simulation experiments (time, parameters, initial value, events)
#' @param p a matrix of parameters with M columns
#' @param abs.tol absolute tolerance, real scalar
#' @param rel.tol relative tolerance, real scalar
#' @param initial.step.size initial value for the step size; the step size will adapt to a value that observes the tolerances, real scalar
#' @param omit an integer that indicates how many of these to omit in this order: fisher information, gradient of the log-likelihood, log-likelihood
#' @return a list of the solution trajectories `y(t;p)` for all experiments (named like the experiments), as well as the output functions
#' @export
#' @keywords ODE
#' @useDynLib uqsa r_gsl_odeiv2_outer_fi
gsl_odeiv2_fi <- function(name,experiments,p,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3, method=0, omit=0){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name) <- so
	}
	if (!file.exists(so)){
            warning(sprintf("[r_gsl_odeiv2_outer] for model name «%s», in directory «%s» file «%s» not found.",name,getwd(),so))
	}
	if (!is.matrix(p)) p <- as.matrix(p)
	y <- .Call(
		r_gsl_odeiv2_outer_fi,
		name,
		experiments,
		p,
		abs.tol,
		rel.tol,
		initial.step.size,
		omit,
		method
	)
	for (i in seq_along(experiments)){
		if ("initialState" %in% names(experiments[[i]]) && length(experiments[[i]]$initialState)==NROW(y[[i]]$state)){
			dimnames(y[[i]]$state) <- list(names(experiments[[i]]$initialState),NULL,NULL)
		}
		if ("outputValues" %in% names(experiments[[i]]) && NCOL(experiments[[i]]$outputValues)==NROW(y[[i]]$func)){
			dimnames(y[[i]]$func) <- list(names(experiments[[i]]$outputValues),NULL,NULL)
		}
	}
	return(y)
}

#' simulates a CRNN ode model with extra work
#'
#' This function calls a C function which solves an initial value
#' problem, derived from a CRNN.
#'
#' @param name the name of the ODE model to simulate (a shared library
#'     of the same name will be dynamically loaded and needs to be
#'     created first)
#' @param experiments a list of `N` simulation experiments (time,
#'     parameters, initial value, events).
#' @param l a matrix of parameters with M columns, in log-space.
#' @param nu a stoichiometry matrix (N×R) where N is the number of state
#'     variables and R the number of reactions, all reactions are
#'     assumed to be reversible.
#' @param m modifiers -- similar to stoichiometry, but indicates whether the
#'     species takes part in the reaction without being consumed.
#' @param abs.tol absolute tolerance, real scalar.
#' @param rel.tol relative tolerance, real scalar.
#' @param initial.step.size initial value for the step size; the step
#'     size will adapt to a value that observes the tolerances, real
#'     scalar.
#' @return a list of the solution trajectories `y(t;p)` for all
#'     experiments (named like the experiments), as well as the output
#'     functions.
#' @export
#' @keywords ODE
#' @useDynLib uqsa r_gsl_odeiv2_outer_CRNN
gsl_odeiv2_CRNN <- function(name,experiments,l,nu,m,abs.tol=1e-6,rel.tol=1e-5,initial.step.size=1e-3,method=0){
	if (is.character(comment(name))){
		so <- comment(name)
	} else {
		so <- paste0("./",name,".so")
		comment(name)<-so
	}
	if (!file.exists(so)){
            warning(sprintf("[r_gsl_odeiv2_outer] for model name «%s», in directory «%s» file «%s» not found.",name,getwd(),so))
	}
	if (!is.matrix(l)) l <- as.matrix(p)
	y <- .Call(r_gsl_odeiv2_outer_CRNN,name,experiments,l,nu,m,abs.tol,rel.tol,initial.step.size,method)
	for (i in seq_along(experiments)){
		if ("initialState" %in% names(experiments[[i]])){
			dimnames(y[[i]]$state) <- list(names(experiments[[i]]$initialState),NULL,NULL)
		}
		if ("outputValues" %in% names(experiments[[i]])){
			dimnames(y[[i]]$func) <- list(names(experiments[[i]]$outputValues),NULL,NULL)
		}
	}
	return(y)
}

#' scrnn returns a closure around gsl_odeiv2_CRNN()
#'
#' the returned value is a function of a variable p that encodes the
#' CRNN in some way. Three user supplied functions are used to extract
#' the three components of a CRNN:
#'
#' - kinetic rate coefficients (in log-space): `l <- parMap(p)`
#' - stoichiometric matrix: `nu <- stoichiometry(p)`
#' - modifier matrix: `m <- modifiers(p)`
#'
#' these three components (one numeric vector, and two matrices) are
#' passed to the simulation procedure. The vector l can be a matrix
#' with M columns. In that case, one simulation per column is
#' performed. The stoichiometry and modifiers remain unchanged
#' throughout.
#' @param experiments list of experiments (inputs are ignored).
#' @param modelName scalar string, can indicate a shared library with
#'     an attached comment attribute.
#' @param parMap (function) extracts kinetic rate coefficients from
#'     its argument.
#' @param stoichiometry (function) extracts the stoichiometry matrix
#'     from its argument.
#' @param modifiers (function) extracts the modifier matrix from its
#'     argument.
#' @param method (integer) integration method key (0:10) corresponds to
#'     these GSL methods: msbdf, msadams, bsimp, rk4imp, rk2imp,
#'     rk1imp, rk8pd, rkck, rkf45, rk4, rk2
#' @return closure that maps one argument (p) to simulation results (y).
#' @export
scrnn <- function(experiments, modelName, parMap=\(p) p$l, stoichiometry=\(p) p$nu, modifiers=\(p) p$m, method = 0){
	N <- length(experiments)
	sim <- function(parMCMC){
		nu <- stoichiometry(parMCMC)
		m <- modifiers(parMCMC)
		l <- parMap(parMCMC)
		yf <- gsl_odeiv2_CRNN(modelName,experiments,l,nu,m,method=0)
		if (N==length(yf)) {
			names(yf) <- names(experiments)
		} else {
			message(sprintf("experiments(%i) should be the same length as simulations(%i), but isn't.",length(experiments),length(yf)))
		}
		return(yf)
	}
	return(sim)
}


#' This creates a closure that simulates the model, similar to simulator.c
#'
#' This is a shorter alternative to simulator.c (C backend). It also
#' returns the log-likelihood, Fisher Information, and the gradient of
#' the log-likelihood, under the assumption that the measurement error
#' is Gaussian.
#'
#' It returns a closure around:
#'     - experiments,
#'     - the model, and
#'     - parameter mapping
#'
#' The returned function depends only on parABC (the sampling
#' parameters). The simulation will be done suing the rgsl backend.
#'
#' This version of the function does not use the parallel package at
#' all and cannot add noise to the simulations.
#'
#' @param experiments a list of experiments to simulate: inital
#'     values, inputs, time vectors, initial times
#' @param modelName a string (with optional comment indicating an .so
#'     file) which points out the model to simulate
#' @param parABC the parameters for the model, subject to change by
#'     parMap.
#' @param parMap the model will be called with parMap(parABC); so any
#'     parameter transformation can happen there.
#' @param method the integration method as an integer (higher numbers are simpler methods, lower numbers are more advanced methods, 0 = msbdf)
#' @param omit integer, omit optional return values, in this order: Fisher Information, gradient of the log-likelihood, the log-likelihood, output functions. Omission includes all previous entries. 'omit = 1' omits only the Fisher Information, omit=3, omits FI, grad-ll, and log-likelihood calculations. 
#' @export
#' @return a closure that returns the model's output for a given
#'     parameter vector, and approximate sensitivity matrices, for
#'     each state variable, function, time-point, and parameter
#'     vector.
#' @examples
#'  #  model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
#'  #  experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
#'  #  parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)
#'
#'  #  modelName <- checkModel("<insert_model_name>_gvf.c")
#'  #  simulate <- simc(experiments, modelName,  parABC)
#'  #  yf <- simulate(parABC)
simfi <- function(experiments, modelName, parMap=identity, method = 0, omit = 0){
	N <- length(experiments)
	## convert data-frames tpo matrices for the C code
	for (i in seq_along(experiments)){
		if (!('data' %in% names(experiments[[i]])) && 'outputValues' %in% names(experiments[[i]])){
			experiments[[i]]$data = t(experiments[[i]][['outputValues']])
		}
		if (!('stdv' %in% names(experiments[[i]])) && 'errorValues' %in% names(experiments[[i]])){
			experiments[[i]]$stdv = t(experiments[[i]][['errorValues']])
		}
	}
	sim <- function(parABC){
		modelPar <- parMap(parABC)
		m <- NCOL(parABC)
		yf <- gsl_odeiv2_fi(
			modelName,
			experiments,
			as.matrix(modelPar),
			method=method,
			omit=min(omit,3)
		)
		if (N==length(yf)) {
			names(yf) <- names(experiments)
		} else {
			message(
				sprintf(
					"experiments(%i) should be the same length as simulations(%i), but isn't.",
					length(experiments),
					length(yf)
				)
			)
		}
		return(yf)
	}
	return(sim)
}

#' This creates a closure that simulates the model
#'
#' This is a shorter alternative to the runModel function (C backend).
#'
#' It returns a closure around:
#'     - experiments,
#'     - the model, and
#'     - parameter mapping
#'
#' The returned function depends only on parABC (the sampling
#' parameters). The simulation will be done suing the rgsl backend.
#'
#' @param experiments a list of experiments to simulate: inital values, inputs, time vectors, initial times
#' @param modelName a string (with optional comment indicating an .so file) which points out the model to simulate
#' @param parABC the parameters for the model, subject to change by parMap.
#' @param parMap the model will be called with parMap(parABC); so any parameter transformation can happen there.
#' @param noise boolean variable. If noise=TRUE, Gaussian noise is added to the output of the simulations. The standard
#'              deviation of the Gaussian noise is equal to the measurement error. If noise=FALSE the output is the
#'              deterministic solution of the ODE system. noise and sensitivity calculations are mutually exclusive.
#' @export
#' @return a closure that returns the model's output for a given parameter vector
#' @examples
#'  #  model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
#'  #  experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
#'  #  parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)
#'
#'  #  modelName <- checkModel("<insert_model_name>_gvf.c")
#'  #  simulate <- simulator.c(experiments, modelName,  parABC)
#'  #  yf <- simulate(parABC)
simulator.c <- function(experiments, modelName, parMap=identity, noise = FALSE, approximateSensitivity = FALSE, method = 0){
	require(rgsl)
	if (noise){
		sim <- function(parABC){
			modelPar <- parMap(parABC)
			yf <- unlist(
				mclapply(
					experiments,
					function(EX) {
						rgsl::r_gsl_odeiv2_outer(modelName, list(EX), as.matrix(modelPar), method=method)
					}
				),
				recursive=FALSE)
			stopifnot(length(experiments)==length(yf))
			for(i in 1:length(experiments)){
				out <- yf[[i]]$func
				l <- dim(out)[2]
				n <- ifelse(is.matrix(parABC),ncol(parABC),1)
				sd <- as.matrix(experiments[[i]]$errorValues)
				if(!is.null(sd)){
					sd[is.na(sd)] <- 0.0
					y <- mclapply(1:n, function(j) {return(out[,,j] + rnorm(l, 0, sd))})
					yf[[i]]$func[1,,] <- do.call(cbind,y)
				}
			}
			return(yf)
		}
	} else if (approximateSensitivity) {
		sim <- function(parABC){
			modelPar <- parMap(parABC)
			yf <- unlist(
				mclapply(
					experiments,
					function(EX) {
						tryCatch(rgsl::r_gsl_odeiv2_outer_sens(modelName, list(EX), as.matrix(modelPar), method=method),
							error = function(e) {prinmt(e); return(NA)})
					}
				), recursive=FALSE)
			stopifnot(length(experiments)==length(yf))
			return(yf)
		}
	} else {
		sim <- function(parABC){
			modelPar <- parMap(parABC)
			yf <- unlist(
				mclapply(
					experiments,
					function(EX) {
						tryCatch(rgsl::r_gsl_odeiv2_outer(modelName, list(EX), as.matrix(modelPar), method=method),
							error = function(e) {print(e); return(NA)})
					}
				), recursive=FALSE)
			stopifnot(length(experiments)==length(yf))
			return(yf)
		}
	}
	return(sim)
}


#' This creates a closure that simulates the model, similar to simulator.c
#'
#' This is a shorter alternative to simulator.c (C backend).
#'
#' It returns a closure around:
#'     - experiments,
#'     - the model, and
#'     - parameter mapping
#'
#' The returned function depends only on parABC (the sampling
#' parameters). The simulation will be done suing the rgsl backend.
#'
#' This version of the function does not use the parallel package at
#' all and cannot add noise to the simulations.
#'
#' @param experiments a list of experiments to simulate: inital
#'     values, inputs, time vectors, initial times
#' @param modelName a string (with optional comment indicating an .so
#'     file) which points out the model to simulate
#' @param parABC the parameters for the model, subject to change by
#'     parMap.
#' @param parMap the model will be called with parMap(parABC); so any
#'     parameter transformation can happen there.
#' @export
#' @return a closure that returns the model's output for a given
#'     parameter vector, and approximate sensitivity matrices, for
#'     each state variable, function, time-point, and parameter
#'     vector.
#' @examples
#'  #  model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
#'  #  experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
#'  #  parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)
#'
#'  #  modelName <- checkModel("<insert_model_name>_gvf.c")
#'  #  simulate <- simc(experiments, modelName,  parABC)
#'  #  yf <- simulate(parABC)
simc <- function(experiments, modelName, parMap=identity, method = 0){
	N <- length(experiments)
	sim <- function(parABC){
		modelPar <- parMap(parABC)
		m <- NCOL(parABC)
		yf <- rgsl::r_gsl_odeiv2_outer_sens(modelName, experiments, as.matrix(modelPar), method=method)
		if (N==length(yf)) {
			names(yf) <- names(experiments)
		} else {
			message(sprintf("experiments(%i) should be the same length as simulations(%i), but isn't.",length(experiments),length(yf)))
		}
		for (i in seq(N)){
			for (j in seq(m)){
				## state variables
				l <- is.finite(yf[[i]]$stateSensitivity[[j]])
				if (any(!l)){
					message(sprintf("state-sensitivity approximation produced %i erroneous elements. Setting invalid elements to 0.0.",sum(!l)))
					##print(yf[[i]]$stateSensitivity[[j]])
					yf[[i]]$stateSensitivity[[j]][!l] <- 0.0
				}
				## functions
				l <- is.finite(yf[[i]]$funcSensitivity[[j]])
				if (any(!l)){
					message(sprintf("function-sensitivity approximation produced %i erroneous elements. Setting invalid elements to 0.",sum(!l)))
					##print(yf[[i]]$funcSensitivity[[j]])
					yf[[i]]$funcSensitivity[[j]][!l] <- 0.0
				}
			}
		}
		return(yf)
	}
	return(sim)
}

#' This creates a closure that simulates the model, similar to simulator.c
#'
#' This is a shorter alternative to simulator.c (C backend).
#'
#' It returns a closure around:
#'     - experiments,
#'     - the model, and
#'     - parameter mapping
#'
#' The returned function depends only on parABC/parMCMC (the sampling
#' parameters). The simulation will be done suing the rgsl backend.
#'
#' This version of the function does not use the parallel package at
#' all and cannot add noise to the simulations. It also doesn't
#' perform sensitivty analysis.
#'
#' @param experiments a list of experiments to simulate: inital
#'     values, inputs, time vectors, initial times
#' @param modelName a string (with optional comment indicating an .so
#'     file) which points out the model to simulate
#' @param parABC the parameters for the model, subject to change by
#'     parMap.
#' @param parMap the model will be called with parMap(parABC); so any
#'     parameter transformation can happen there.
#' @export
#' @return a closure that returns the model's output for a given
#'     parameter vector
#' @examples
#'  #  model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
#'  #  experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
#'  #  parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)
#'
#'  #  modelName <- checkModel("<insert_model_name>_gvf.c")
#'  #  simulate <- simcf(experiments, modelName,  parABC)
#'  #  yf <- simulate(parABC)
simcf <- function(experiments, modelName, parMap=identity, method = 0){
	N <- length(experiments)
	sim <- function(parABC){
		modelPar <- parMap(parABC)
		m <- NCOL(parABC)
		yf <- rgsl::gsl_odeiv2_outer(name,experiments,modelPar,abs.tol,rel.tol,initial.step.size, method)
		return(y)
	}
	return(sim)
}

#' checkModel tries to establish the simulation file for a given model
#'
#' This function returns the model name, with some additional comments
#' about the file
#'
#' As an alternative to this function, it is sufficient to write
#'
#' ```
#' modelName <- "test_ode_model"             # or some other model name
#' comment(modelName) <- "test_ode_model.so"
#' ```
#'
#' This function will not attempt to find a model file, other than in
#' the current directory. But, `checkModel` will compile a GSL
#' compatible C source file into a shared object if `modelFile` ends
#' with `.c` and stop if that doesn't work. The compiler is called
#' using a system call, which may be incorrect for your system -- if
#' this funciton fails, you'll have to make the shared library of the
#' model using the correct compiler and options for your system.
#'
#' In any case, this function stops execution if the model file
#' doesn't exist, because simulations are not possible.
#'
#' @export
#' @param modelName a string
#' @param modelFile a string, if the model file is different from
#'     "modelName.R". If the file name ends in .c, the c source will be
#'     compiled to a shared library.
#' @return modelName with an additional comment about which file to use for simulations
checkModel <- function(modelName,modelFile=NULL){
	if (is.null(modelFile)) {
		modelFile <- paste0(modelName,c('.R','_gvf.c','.so'));
		modelFile <- modelFile[file.exists(modelFile)]
		stopifnot(!is.null(modelFile) && length(modelFile)>0)
		modelFile <- modelFile[1]
	}
	if (grepl('[.]c$',modelFile,useBytes=TRUE)){
		stopifnot(file.exists(modelFile))
		message('building a shared library from c source, and using GSL odeiv2 as backend (pkg-config is used here).')
		LIBS <- "`pkg-config --libs gsl`"
		CFLAGS <- "-shared -fPIC `pkg-config --cflags gsl`"
		so <- sprintf("./%s.so",modelName)
		command_args <- sprintf("%s -o '%s' '%s' %s",CFLAGS,so,modelFile,LIBS)
		message(paste("cc",command_args))
		system2("cc",command_args)
		stopifnot(file.exists(so))
		comment(modelName)<-so
	} else if (grepl('[.]so$',modelFile,useBytes=TRUE)) {
		stopifnot(file.exists(modelFile))
		message(sprintf('Will use pre-existing %s for simulations.',modelFile))
		comment(modelName) <- modelFile
	} else if (grepl('[.]R$',modelFile,useBytes=TRUE)){
		message(sprintf('will use the %s for simulations (deSolve backend)',modelFile))
		comment(modelName) <- modelFile
	}
	return(modelName)
}

#' default distance function for one experiment
#'
#' if each experiment corresponds to one simulation and is fully
#' quanitified by itself, then calculating the overall distance
#' between data and experiment can be done one by one. This function
#' describes the default way a simulation is compared to data.
#'
#' If the data is more complex, and two or more simulations are needed
#' to calculate one distance value then the objective-Function needs
#' to be entirely user-supplied. This is the case with experiments
#' that have a "control" -- this is needed when the measurement is in
#' arbitrary units and only makes sense comparatively to a secondary
#' (control) scenario.
#'
#' This function will be used if none is provided by the user.
#'
#' The funcSim values need to be supplied as a matrix of size N×T with
#' N the length of the model's output vectors and T the amount of
#' measurement times (this is how the rgsl package returns the
#' simulation results).
#' @param funcSim a matrix, contains model solution (output values),
#'     columns of output vectors
#' @param dataVAL a data.frame of experimental data
#' @param dataERR a data.frame of measurement errors, if available,
#'     defaults to the maximum data value.
defaultDistance <- function(funcSim,dataVAL,dataERR=max(dataVAL)){
	if (all(is.finite(funcSim))){
		distance <- mean(abs(funcSim-t(dataVAL))/t(dataERR), na.rm=TRUE)
	} else {
		distance <- Inf
	}
	return(distance)
}

#' default ABC acceptance probability function for one experiment
#'
#' if each experiment corresponds to one simulation and is fully
#' quanitified by itself, then calculating the overall distance
#' between data and experiment can be done one by one. This function
#' describes the default way a simulation is compared to data.
#'
#' If the data is more complex, and two or more simulations are needed
#' to calculate one distance value then the objective-Function needs
#' to be entirely user-supplied. This is the case with experiments
#' that have a "control" -- this is needed when the measurement is in
#' arbitrary units and only makes sense comparatively to a secondary
#' (control) scenario.
#'
#' This function will be used if none is provided by the user.
#'
#' The funcSim values need to be supplied as a matrix of size N×T with
#' N the length of the model's output vectors and T the amount of
#' measurement times (this is how the rgsl package returns the
#' simulation results).
#' @param funcSim a matrix, contains model solution (output values),
#'     columns of output vectors
#' @param dataVAL a data.frame of experimental data
#' @param dataERR a data.frame of measurement errors, if available,
#'     defaults to the maximum data value.
defaultAcceptance <- function(funcSim,dataVAL,dataERR=max(dataVAL)){
	n <- prod(dim(funcSim))
	if (all(is.finite(funcSim))){
		ABCP <- exp(-0.5*sum((abs(funcSim-t(dataVAL))/t(dataERR))^2, na.rm=TRUE))#/(sqrt(2*pi)^n*prod(as.numeric(dataERR)))
	} else {
		ABCP <- 0.0
	}
	return(ABCP)
}


#' creates Objective functions from ingredients
#'
#' the returned objective function has only one argument: the ABC
#' variables that shall be mapped to ODE-model parameters.
#'
#' @export
#' @param experiments a list of simulation experiments
#' @param modelName and model storage file as comment
#' @param distance a function that calculates ABC scores
#' @param parMap a function that transforms ABC variables into acceptable model parameters
#' @param simulate closure that simulates the model
#' @return an objective function
makeObjective <- function(experiments,modelName=NULL,distance,parMap=identity,simulate=NULL)
{
	if (is.null(simulate) && !is.null(modelName)){
		simulate <- function(par){
			return(runModel(experiments, modelName,  par, parMap))
		}
	}
	N <- length(experiments)
	Objective <- function(parABC){
		out <- simulate(parABC)
		n <- ifelse(is.matrix(parABC),ncol(parABC),1)
		S <- matrix(Inf,nrow=N,ncol=n)
		rownames(S) <- names(experiments)
		if (is.null(out)) return(S)
		for(i in 1:N){
			if (!is.null(experiments[[i]]) && !is.null(out[[i]])){
				S[i,] <- unlist(mclapply(1:n, function(j) distance(out[[i]]$func[,,j], experiments[[i]]$outputValues, experiments[[i]]$errorValues)))
			}
		}
		return(S)
	}
	return(Objective)
}
