# uqsa 0.8.1

- Fixed a bug in `parse_concise()`.
- Reduced the amount of warnings about _possible_ model feature
  incompatibility discovered while parsing the model.

In the future, this feature analysis will issue warnings only when
it is a true conflict, not a possible conflict (e.g.: scheduled events
modifying species that are subject to conservation law experssions)
