# Updates the named values of vector `v` with values mentioned in data.frame d

Given a named vector, this function will create a mtarix, with several
copies of the vector, with different values, dependeing on some context.
The context is each row of a `data.frame`, with columns referring to the
names of `v`, and values for the entries in `v` for that row's context.
Example: the rows can be different experinemnts, with each experiment
assigning new values to some members of `v` (but not necessarily all).

## Usage

``` r
update_values(v, d, as_type = "numeric")
```

## Arguments

- v:

  a named vector

- d:

  data.frame with column names that correspond to those of `v`

- as_type:

  a character scalar indicating a type
  ('character','numeric','logical',etc.)

## Value

a matrix of dimension: length(v) × NROW(d)

## Examples

``` r
if (FALSE) { # \dontrun{
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
iv <- values(m$Compound)
IV <- update_values(iv,m$Experiments)
print(IV)
} # }
```
