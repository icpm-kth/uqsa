# Finds the values in a data.frame

given a data.frame, this function will look for a column that contains
some kind of value and retrieve it. The returned numeric vector will be
named.

## Usage

``` r
values(df)
```

## Arguments

- df:

  a data frame with a "value" column

## Value

a named numeric vector

## Details

If the values contain a standard error, the returned value is of class
*errors*.
