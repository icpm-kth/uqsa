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
#'     e.g. "3 cm"
#' @param target string, target unit, e.g. "m"
#' @return a numeric value y: val*originalUnit = y*targetUnit, the
#'     target unit is attached to the returned value, as a comment.
#' @examples
#' y <- "21 cm" %as% "inches"
`%as%` <- function(txtUnit,target){
	if (system2("command",args=c("-v","units"),stderr=FALSE,stdout=FALSE)==0){
		f <- as.double(
			system2(
				"units",
				args=c(
					paste0("--",c("strict","compact")),
					"-1",
					sprintf("'%s'",txtUnit),
					sprintf("'%s'",target)
				),
				stdout=TRUE
			)
		)
		comment(f) <- target
	} else {
		warning("The 'units' utility must be installed (system program, not R).")
	}
	return(f)
}
