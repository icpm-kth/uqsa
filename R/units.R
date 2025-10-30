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
