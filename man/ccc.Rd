% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gillespie.R
\name{ccc}
\alias{ccc}
\title{Concentration to Count Conversion}
\usage{
ccc(sb)
}
\arguments{
\item{sb}{SBtab list of data.frames with one of the members called
'Compound', with a mandatory "!Unit" column.}

\item{LV}{Avogadro's constant multiplied by system volume}
}
\value{
f a conversion factor: n <- f*x, where n is a particle
count and x a concentration in the specified volume.
}
\description{
This function returns the factor f for converting a concentration
to particle numbers (counts) based on the units specified in the
SBtab file and the volume of the system encoded as LV (Avogadro's
constant × Volume).
}
\details{
The implicit assumption is that the LV constant the user passes to
this function has a compatible unit of volume to the SBtab unit of
concentration. LV should be in litres and the concentrations
also in litres up to SI prefixes. All of these are OK: M, nmol/l,
mol/dl, kilomol/megalitre.

If the concentrations are 'something per cubic metre', then LV has to be in m³
as well. The unit of LV is not checked in any way. But the unit of
concentration will be parsed an converted to mol/l.
}
