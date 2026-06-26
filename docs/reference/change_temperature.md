# Should 2 Markov chains exchange their temperatures

This function makes a Boolean choice about chnages in temperature, based
on the log(liklihood) values of two Markov chains in a parallel
tempering setting.

## Usage

``` r
change_temperature(b1, ll1, b2, ll2)
```

## Arguments

- b1:

  the inverse temperature of chain 1

- ll1:

  the log-likelihood of chain 1

- b2:

  the inverse temperature of chain 2

- ll2:

  the log-likelihood of chain 2

## Value

TRUE is the chains should swap their temperatures

## Details

This function is useful if `mpi.send()` and `mpi.recv()` are used.

## Examples

``` r
b <- c(1.0,0.5)
if (change_temperature(b[1],-850,b[2],-600)){ # with some randomness
  message(sprintf("yes, swapping temperature %f <=> %f",b[1],b[2]))
} else {
  message(sprintf("no, temperature %f, and %f stay unchanged",b[1],b[2]))
}
#> yes, swapping temperature 1.000000 <=> 0.500000
```
