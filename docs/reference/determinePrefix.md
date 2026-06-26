# Determine a prefix from a character vector str of similar contents

The result is such that `all(startsWith(str,determinePrefix(str)))` is
`TRUE`.

## Usage

``` r
determinePrefix(str, split = "-", collapse = "-")
```

## Arguments

- str:

  a character vector

- split:

  the token to use for strsplit instead of '-', this should be
  `character(0)` if you want to split letter by letter

- collapse:

  the words constituents in the input that are found to be uniform in
  the input are connected via [paste](https://rdrr.io/r/base/paste.html)
  and this "collapse" value.

## Value

the prefix common to all entries of str.

## Details

By default, the strings are assumed to be '-' separated words, and a
series of words is found to be the prefix if all entries start with that
set of words.

The normal case is c("abc-1","abc-2b","abc-2a") maps to "abc"

## Examples

``` r
files <- sprintf("smmala-sample-%i-of-3.RDS",seq(1,3))
pref <- determinePrefix(files)
print(pref)
#> [1] "smmala-sample"
```
