# Changelog

## uqsa 0.8.1

CRAN release: 2026-10-02

- Fixed a bug in
  [`parse_concise()`](https://icpm-kth.github.io/uqsa/reference/parse_concise.md).
- Reduced the amount of warnings about *possible* model feature
  incompatibility discovered while parsing the model.

In the future, this feature analysis will issue warnings only when it is
a true conflict, not a possible conflict (e.g.: scheduled events
modifying species that are subject to conservation law experssions)
