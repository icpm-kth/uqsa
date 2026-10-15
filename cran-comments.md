## uqsa 0.8.3

- The ODE solver was made more robust with respect to missing data
  values, and _standard errors_ that are specified as `NA`, or `NaN`,
  or `Inf` (because of the missing data values).
- Some functions were declared _internal_ via keyword, to not clutter
  the packages Index on the `help(package='uqsa')` page
- Some function exports were removed (hard internal), for functions
  that were of dubious (but possible) usefulness to the user
- Made the title of each manual page more similar, usually starting
  with a verb now (adjusted to work well with the function Index)
- the `%as%` function was made to work with both _GNU units_ and
  _FreeBSD units_, which is the default on MACOS, it now uses only
  options available in both variants

Eventually the `%as%` function will be re-implemented to either use a
library, or custom code. The tool use is an interim solution.

