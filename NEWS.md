# uqsa 0.8.4

- added a print function for the model itself
- the model printout is organized in box-art to condense the information

# uqsa 0.8.3

- The ODE solver was made more robust with respect to missing data
  values, and _standard errors_ that are specified as `NA`, or `NaN`,
  or `Inf` (because of the missing data values)
- Some functions were declared _internal_ via keyword, to not clutter
  the packages Index on the `help(package='uqsa')` page
- Some function exports were removed, for functions
  that were of dubious (but possible) usefulness to users
- Made the title of each manual page more similar, usually starting
  with a verb now (adjusted to work well with the function Index)
- the `%as%` function was made to work with both _GNU units_ and
  _FreeBSD units_, which is the default on MACOS, it now uses only
  options available in both variants


# uqsa 0.8.1

- Fixed a bug in `parse_concise()`.
- Reduced the amount of warnings about _possible_ model feature
  incompatibility discovered while parsing the model.

In the future, this feature analysis will issue warnings only when
it is a true conflict, not a possible conflict (e.g.: scheduled events
modifying species that are subject to conservation law experssions)
