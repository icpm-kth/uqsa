# Creates a Matrix with reaction kinetic entries

Given a data.frame that describes reaction, this function extracts the
forward flux and backward flux from the kinetic.law column.

## Usage

``` r
flux_matrix(Reaction)
```

## Arguments

- Reaction:

  a data.frame with a column named 'kinetic.law'
