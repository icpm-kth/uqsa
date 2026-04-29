#' Construct Code
#'
#' Interpret the first argument and generate code in the specified
#' language for the model type.
#'
#' Whenever the model `Model` is of type `"cme"`, the LV parameter is
#' used to determine the actual number of molecules in the
#' system. Otherwise it is ignored.
#'
#' @param Model either CME or ODE model
#' @param language either C or R
#' @param LV Avogadro's constant multiplied by the system's volume in
#'     litres, only used for CME models
#' @export
#' @return a character vector with the code
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' C <- generate_code(o)
#' cat(head(C),sep="\n")
generate_code <- function(Model,language="C", LV=6.02214076e+8){
	if (is(Model,"ode") && tolower(language)=="c"){
		C <- generateCode(Model)
	} else if (is(Model,"ode") && tolower(language)=="r"){
		C <- generateRCode(Model)
	} else if(is(Model,"cme") && tolower(language)=="c"){
		C <- generateGillespieCode(Model,LV)
	} else {
		stop(sprintf("Requested source-code type not implemented (yet): %s in %s",class(Model),as.character(language)))
	}
	comment(C) <- Model$name
	return(C)
}

#' Compile C code to shared library
#'
#' Calls `R CMD SHLIB` to create the model's shared library.
#'
#' The first argument can be a raw character scalar with just the path
#' of the c code to be compiled, or alternatively an object that has
#' this information stored within it. The models returned by [as_cme]
#' and [as_ode] can both carry this information, attach it via
#' [c_path<-].
#'
#' @param file the c file that is to be compiled, OR an ODE/CME object
#'     with a c.file defined and recorded in it.
#' @export
#' @return the path of the created shared library
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' C <- generate_code(o)
#' c_path(o) <- write_c_code(C)
#' so_path(o) <- shlib(o)
#' print(o)
#' if (file.exists(so_path(o))) cat("shared library exists.\n")
shlib <- function(file){
	if (is(file,"ode")) {
		so_name <- file$name
		file <- c_path(file)
	} else if (is(file,"cme")) {
		so_name <- file$name
		file <- c_path(file)
	} else {
		so_name <- sub("[.]c$","",basename(file))
	}
	if (is.null(file)) stop("No C file specified")
	if (!file.exists(file)) stop(sprintf("%s does not exist.",file))

	so <- file.path(dirname(file), paste0(so_name, .Platform$dynlib.ext))
	cflags <- system2("pkg-config", c("--cflags", "gsl"), stdout = TRUE)
	libs <- system2("pkg-config", c("--libs", "gsl"), stdout = TRUE)
	compile_env <- c(
		if (nzchar(cflags)) sprintf("PKG_CPPFLAGS='%s'",cflags),
		if (nzchar(libs)) sprintf("PKG_LIBS='%s'", libs)
	)
	status <- system2(
		command = file.path(R.home("bin"), "R"),
		args = c("CMD", "SHLIB",file),
		env = compile_env,
		stdout = TRUE,
		stderr = TRUE,
		wait = TRUE
	)
	if (!file.exists(so)) {
		cat(status)
		warning(sprintf("Building %s failed.",so))
	}
	return(so)
}



#' Write the C code to a file
#'
#' This function does not compile the code, it only writes it to a
#' file in a temporary location (tempdir). By default, the name of the
#' file will contain the hash of the entire code.
#'
#' If instead of a character vector, an ode or cme object is passed,
#' this function will generate code from it with default options.
#'
#' @param C the code to write, as a character array.
#' @param model.name a string with no special characters, will be used in the file name
#' @param file override the default file name (based on hashing)
#' @export
#' @return the path of the written file
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' C <- generate_code(o)
#' c_path(o) <- write_c_code(C)
#' print(o)
#' if (file.exists(c_path(o))) cat("c file exists.\n")
write_c_code <- function(C, model.name=comment(C), file=file.path(tempdir(),digest::digest(C,"xxhash64"),paste0(model.name,".c"))){
	if (is(C,"ode")) { # an ode model was passed instead of code
		ode <- C
		model.name <- ode$name
		C <- generate_code(ode)
	} else if (is(C,"cme")){
		cme <- C
		model.name <- cme$name
		C <- generate_code(cme)
	}
	cat(sprintf("Writing file: %s\n",file))
	if (!dir.exists(dirname(file))){
		dir.create(dirname(file),recursive=TRUE)
	}
	cat(C,sep="\n",file=file)
	return(file)
}

#' Writes code to file and compiles
#'
#' This function accepts an ode model, or cme model, generates code,
#' compiles it to a shared library, and returns a changed object.
#' possibly changed by the user. It writes the contents to a c file
#' named 'modelName_gvf.c'. This file is compiled to './modelName.so'
#' using normal command line tools, not `R CMD SHLIB`
#'
#' This entire function can be replaced with a call to `cat()` and
#' then compiling the written file in the system's shell.
#'
#' @export
#' @param M ode or cme Model for which code is generated and written to a file
#' @return a copy of `o` with file paths added to it
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' print(o)
write_and_compile <- function(M){
	C <- generate_code(M)
	f <- tempfile(pattern=sprintf("%s_",M$name), fileext=".c")
	cat(C,sep="\n",file=f)
	c_path(M) <- f
	so_path(M) <- shlib(f)
	return(M)
}
