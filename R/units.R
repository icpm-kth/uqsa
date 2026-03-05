#' find unit category
#'
#' This function reads a unit, without SI prefix, and returns a string
#' from a smaller subset of unit kinds, similar to what is defined in
#' SBML. This normalizes the various ways to write the same unit:
#' "meter", "m" and "metre".
#'
#' The unit kind of "m" is "metre", the kind of "g" is "gram".
#'
#' @param kind the unnormalized string that humans use to write a unit
#'     (but without prefix)
#' @return normalized category name of the unit kind: litre, mole,
#'     metre, kilogram, gram, ampere, candela, second, kelvin, hour,
#'     molarity, dimensionless. defaults to dimensionless.
#' @examples
#' >  unit.kind("meter")
#' [1] "metre"
unit.kind <- function(kind){
	stopifnot(is.character(kind) && length(kind)==1)
	if (grepl("^(l|L|litre|liter)$",kind)){
		k <- "litre"
	} else if (grepl("^(mole?)$",kind)) {
		k <- "mole"
	} else if (grepl("^(m|meter|metre)$",kind)) {
		k <- "metre"
	} else if (grepl("^(kg|kilogram)$",kind)){
		k <- "kilogram"
	} else if (grepl("^(g|gram)$",kind)){
		k <- "gram"
	} else if (grepl("^(A|ampere)$",kind)){
		k <- "ampere"
	} else if (grepl("^(cd|candela)$",kind)){
		k <- "candela"
	} else if (grepl("^(s|second)$",kind)){
		k <- "second"
	} else if (grepl("^(K|kelvin)$",kind)){
		k <- "kelvin"
	} else if (grepl("^(N|[Nn]ewton)$",kind)){
		k <- "newton"
	} else if (grepl("^(h|hour)$",kind)){
		k <- "hour"
	} else if (grepl("^(M|molarity)$",kind)){
		k <- "molarity"
	} else {
		stop(sprintf("The unit kind «%s» is not known, yet."))
	}
	return(k)
}

#' Unit scale from SI prefix
#'
#' This function reads a prefix from a string and returns the exponent
#' (base-10) that this prefix represents.
#'
#' @param prefix a string, e.g.: "M", "mega", "m", "milli", "µ", "micro", etc.
#' @return an integer that corresponds to the prefix, defaults to 0.
#' @examples
#' > unit.scale("M")
#' [1] 6
#' > unit.scale("µ")
#' [1] -6
unit.scale <- function(prefix){
	stopifnot(is.character(prefix) && length(prefix)==1)
	if (grepl("^G$|^giga$",prefix)){
		s <- 9
	} else if (grepl("^M$|^mega$",prefix)){
		s <- 6
	} else if (grepl("^k$|^kilo$",prefix)){
		s <- 3
	} else if (grepl("^h$|^hecto$",prefix)){
		s <- 2
	} else if (grepl("^d$|^deci$",prefix)){
		s <- -1
	} else if (grepl("^c$|^centi$",prefix)){
		s <- -2
	} else if (grepl("^m$|^milli$",prefix)){
		s <- -3
	} else if (grepl("^u$|^µ$|^\xCE\xBC$|^micro$",prefix)){
		s <- -6
	} else if (grepl("^n$|^nano$",prefix)){
		s <- -9
	} else if (grepl("^p$|^pico$",prefix)){
		s <- -12
	} else if (grepl("^f$|^femto$",prefix)){
		s <- -15
	} else {
		s <- 0
	}
	return(s)
}

#' Converts a unit to a string that works as an identifier
#'
#' Some formats require a name for a unit definition. This functions
#' creates a name from a unit, converting math/symbols to text. The
#' returned value should work as an SBML unit id.
#'
#' @export
#' @param unit.str the original string representastion of that unit
#' @param prnt logical switch: if TRUE, the name will be printed.
#' @return unit.id string
#' @examples
#' > unit.id("s^9")
#' [1] "s_to_the_power_of_9"
#'
#' > unit.id("cm^2")
#' [1] "cm_square"
#'
#' > unit.id("1/s")
#" [1] "one_over_s"
unit.id <- function(unit.str,prnt=FALSE){
	uid <- unit.str
	uid <- sub("^1$","dimensionless",uid)
	uid <- gsub("1/","one_over_",uid)
	uid <- gsub("/","_per_",uid)
	uid <- gsub("[*[:blank:]]","_",uid)
	uid <- gsub("[()]","",uid)
	uid <- gsub("\\^2","_square",uid)
	uid <- gsub("\\^3","_cube",uid)
	uid <- gsub("\\^([0-9]+)","_to_the_power_of_\\1",uid)
	uid <- make.names(uid,unique=FALSE)
	if (prnt){
		message("units in «!Unit» column:")
		print(unit.str)
		message("automatically created sbml unit ids:")
		print(uid)
	}
	return(uid)
}

`%~%` <- function(text, pattern){
	r <- regexec(pattern=pattern,text=text)
	m <- regmatches(text,r)
	return(m)
}


#' Simple unit from string
#'
#' This function takes a simple, human readable unit (without '*' or '/'), from a string
#' and returns a data.frame with the unit's meaning.
#'
#' In this context, a simple unit is just a prefix, a unit kind, and an exponent, e.g. cm^2
#' A not-simple unit is: m/s, kg*m/s^2, kg*h
#' @param u a unit with no fractions or products
#' @return a data.frame with the unit's properties
simple.unit <- function(u=NULL){
	## defaults
	u.m <- 1
	u.x <- 1
	u.s <- 0
	u.k <- "dimensionless"
	## an empty unit means that the value is dimensionless (the unit is '1')
	if (!nzchar(u)) return(data.frame(scale=u.s,multiplier=u.m,exponent=u.x,kind=u.k))
	## um, actually, kg is an SI unit "kind", but doesn't take other prefixes
	prefix.pattern <- "(G|giga|M|mega|k|kilo|h|hecto|c|centi|m|milli|u|\xCE\xBC|\xc2\xb5|micro|n|nano|p|pico|f|femto)?"
	unit.name.pattern <- "(l|L|liter|litre|g|gram|mole?|h|hour|s|second|m|meter|metre|K|kelvin|cd|candela|A|ampere|M|molarity|N|[Nn]ewton)"
	exponent.pattern <- "\\^?([-+]?[0-9]+)?"
	pat <- paste0("^",prefix.pattern,unit.name.pattern,exponent.pattern,"$")
	if (grepl("^kg|kilogram$",u)){
		u.k <- "kilogram"
		u.s <- 0
		u.x <- 1
	} else {
		m <- unlist(u %~% pat)
		if (length(m) > 0){
			u.s <- unit.scale(m[2])
			u.k <- unit.kind(m[3])
			if (nchar(m[4])>0) u.x <- as.numeric(m[4])
		}
	}
	## some special units that need fixing
	if (u.k == "hour") {
		u.k  <- unit.kind("s")
		u.m <- 60
	} else if (u.k == "molarity") {
		u.k <- c(unit.kind("mole"),unit.kind("litre"))
		u.m <- c(u.m,1)
		u.x <- c(u.x,-u.x)
		u.s <- c(u.s,0)
	}
	return(data.frame(scale=u.s,multiplier=u.m,exponent=u.x,kind=u.k))
}

trimmed_split <- function(a,b,fixed=TRUE,...){
	s <- trimws(unlist(strsplit(a,"/",fixed=fixed,...)))
	l <- nzchar(s)
	return(s[l])
}

#' Unit Interpreter
#'
#' This function will try its best to interpret strings like
#' "liter/(nmol ms)"
#' rules: 1. only one slash is allowed
#'        2. M can be mega or mol/l: writing M for molarity will treat
#'           molarity as it's own unit kind; writing "molarity"
#'           will be translated into two SI units (mol and litre)
#'        3. prefixes and units can be words or single letters
#'        4. everything after a slash is the denominator
#'        5. u is an accepted replacement for μ
#'           (unicode greek mu or unicode micro symbol)
#'        6. no parentheses (ignored): "(m/s)*kg" will be misinterpreted
#'
#' this retruns a data.frame with components as in the sbml standard:
#' kind, multiplier, scale and exponent since there is only one
#' slash,parentheses do nothing everything after a slash is the
#' denominator, so: l/mol s is the same as (l)/(mol s) Remark: not all
#' units are understood.
#' @param unit.str a string that contains a human readable unit
#' @param verbose if TRUE, this function prints what it does (to find
#'     problems)
#' @return data.frame with an interpretation of the unit (multiplier
#'     is unused here, but may be used later to deal with units such
#'     as hours (kind=second, multiplier=60)
#' @export
#' @examples
#' > unit.from.string("m/s")
#'   scale multiplier exponent   kind
#' 1     0          1        1  metre
#' 2     0          1       -1 second
#'
#' > unit.from.string("micromolarity")
#'   scale multiplier exponent  kind
#' 1    -6          1        1  mole
#' 2     0          1       -1 litre
#'
#' > unit.from.string("µM")
#'   scale multiplier exponent     kind
#' 1    -6          1        1 molarity
unit.from.string <- function(unit.str){
	if (!is.character(unit.str)){
		print(unit.str)
		stop("unit.str has to be a charcter vector of length 1.")
	}
	stopifnot(length(unit.str)==1)
	unit <- NULL
	a <- gsub("[()]","",unit.str)
	a <- gsub("molarity","mol l^-1",a);
	if (grepl("/",unit.str)){
		a <- trimmed_split(a,"/")
	}
	n <- length(a)
	stopifnot(n==1 || n==2)
	for (j in 1:n){
		b <- trimmed_split(a[j],"[* ]",fixed=FALSE)
		for (u in b){
			su <- simple.unit(u)
			if (j>1) su$exponent <- -su$exponent
			unit <- rbind(unit,su)
		}
	}
	comment(unit) <- unit.str
	attr(unit,'id') <- unit.id(unit.str)
	return(unit)
}

#' converts a unit data.frame into a printable string
#'
#' This is a crude function to make a printable representation of a
#' unit data.frame
#'
#' This function is similar to unit.info.
#' @export
#' @param unit a data.frame created by unit.from.string()
#' @return a string representation of that data.frame purely for printing
unit_as_character <- function(unit){
	if (any(unit$multiplier!=1.0)){
		return(
			paste(
				sprintf(
					"(%g * %s * 10^(%i))^(%i)",
					unit$multiplier,
					unit$kind,
					unit$scale,
					unit$exponent
				),
				collapse='*'
			)
		)
	} else if (any(unit$exponent!=0)){
		return(
			paste(
				sprintf(
					"(%s*10^(%i))^(%i)",
					unit$kind,
					unit$scale,
					unit$exponent
				),
				collapse='*'
			)
		)
	} else {
		return(
			paste(
				sprintf(
					"(%s*10^(%i))",
					unit$kind,
					unit$scale
				),
				collapse='*'
			)
		)
	}
}

#' Prints an interpretation string of a unit
#'
#' given a string describing a unit of measurement, this function
#' prints the interpretation on screen, rather than returning it as a
#' data.frame
#'
#' @param unit.str unit string
#' @param unit optionally, the data.frame that describes the unit
#' @export
#' @examples
#' > unit.info("km/h",unit.from.string("km/h"))
#' «km/h» has been interpreted as:
#'	(1 × metre×10^(3))^(1)
#'	(60 × second×10^(0))^(-1)
unit.info <- function(unit.str,unit=unit.from.string(unit.str)){
	Info <- sprintf(
		"(%g × %s×10^(%i))^(%i)",
		unit$multiplier,
		unit$kind,
		unit$scale,
		unit$exponent
	)
	cat(
		sprintf("«%s» has been interpreted as: ",unit.str),
		Info,
		sep='\n'
	)
	return(Info)
}

#' %as% is a binary operator on strings with units in them
#'
#' The function calls the units utility and converts the string on the
#' left into the unit on the right, e.g.: "cm" %as% "inches", both
#' units can contain numbers. Any input that is accepted by the units
#' utility is acceptable, as long as it makes sense with the command
#' line arguments: `units --strict --compact -1 "$originalUnit" "$targetUnit"`
#'
#' @export
#' @param txtUnit a string with numeric values, including units,
#'     e.g. "3 cm", can be a character vector
#' @param target string, target unit, e.g. "m", must be scalar
#' @return a numeric value y: val*originalUnit = y*targetUnit, the
#'     target unit is attached to the returned value, as a comment.
#' @examples
#' y <- "21 cm" %as% "inches"
#' y <- "12 μmol/L" %as% "mol/L"
#' print(comment(y))
#' y <- "12 mol/m^3" %as% "mmol/L"
`%as%` <- function(txtUnit,target){
	if (system2("command",args=c("-v","units"),stderr=FALSE,stdout=FALSE)){
		warning("The 'units' utility must be installed (system program, not R).")
		return(NA)
	}
	if (length(target)!=1) warning("There must be exactly one target unit.")
	f <- as.numeric(
		sapply(
			txtUnit,
			\(u) return(
				system2(
					command="units",
					args=c(
						paste0("--",c("strict","compact")),
						"-1",
						sprintf("'%s'",as.character(u)),
						sprintf("'%s'",as.character(paste(target,collapse="")))
					),
					stdout=TRUE
				)
			)
		)
	)
	attr(f,"unit") <- target
	names(f) <- txtUnit
	return(f)
}
