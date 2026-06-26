# Reverse look-up of method name from key

These are the methods in the gsl library (documented in the official
documentation), but in reverse order, as they are approximately ordered
by complexity, with more complex methods usually being better (but
slower).

## Usage

``` r
name_method(key = seq(0, 10))
```

## Arguments

- key:

  an integer from 0 to 10 (this is used as an offset in c, for 11 items)

## Value

a string representation of the integration method.

## Details

It is therefore a reasonable approach to try methods from the more
complex end of the list first and try the next method if the solutions
are too slow. But we need to check the accuracy/stability of the result.
The mapping between method names and keys:

         msbdf:  0
       msadams:  1
         bsimp:  2
        rk4imp:  3
        rk2imp:  4
        rk1imp:  5
         rk8pd:  6
          rkck:  7
         rkf45:  8
           rk4:  9
           rk2: 10

The returned value is an integer index.

## Examples

``` r
print(name_method())
#>  [1] "msbdf"   "msadams" "bsimp"   "rk4imp"  "rk2imp"  "rk1imp"  "rk8pd"  
#>  [8] "rkck"    "rkf45"   "rk4"     "rk2"    
```
