## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

### Changes since last manual review

Most recent (and first) human review was by Konstanze Lauseker <konstanze.lauseker@wu.ac.at> on the 7th of June 2026. Thank you!

Summary of what we did to address the points raised by the review:

- moved all S3 object print functions (`print.*`) to their own file:
  R/print.R
- all other functions have guarded print/message/cat statements that
  only trigger when a `verbose` option is set, or if it is part of a
  `warning` or `stop` (printing an R object)
- converted most print and cat messages to actual `message` calls,
  unless they print a complex R object (e.g `data.frame`)
- every example, vignette, and article that changes user settings
  restores the original `options` and `par` settings
- checked that `on.exit()` is used where appropriate (in functions)
- replaced all acronyms with the fully spelled out meaning, as we use
  each acronym only once in the description
- all functions have a `@return` (`\value`) field
- all functions with examples are either exported or have no Rd file
- the only remaining function that has a `\dontrun` example is
  `mcmc_mpi`, it cannot run outside of an MPI context
   + started with `mpirun`
   + uses 'pbdMPI' for MPI interoperation
- removed all `donttest` occurrences

### Detailed Explanations of Changes

We moved all `print.*` methods to their own file `R/print.R`, to make
it _easier_ to check for stray cats and prints.

Some examples used to be guarded by `\dontrun` _only_ because it was
difficult to make them work in under 5 seconds. We replaced all of
these with an `if (interactive()) {...} else {...}` guard for the most
work-intensive function call. The setup for the example is the same in
either if-case.  When run interactively, the examples will use a slightly
_bigger_ problem size and run for 10s to 20s (approximately). When run
in a non-interactive session, the example is reduced to the absolute
_smallest_ problem size to stay under 5s.

All examples that would produce plots are also guarded by
`interactive()`.

All Rd files have a `\value`, we check using: `grep -c '\\value' man/*.Rd`

We check for stray `print()` calls using: `grep '^[[:alnum:][:blank:]]*print(' *.R`
All print calls that we could find have a justification comment, like this:
- `print(...) # guarded by verbose`
- `print(...) # part of warning`

Similarly for stray cats: `grep '^[[:alnum:][:blank:]]*cat(' *.R`. Remaining occurrences are commented, e.g.:
- `cat( # writes to a file`
