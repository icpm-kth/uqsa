# Returns information about parameter conversion

Given parameters of a reaction kinetic coefficient model with
concentrations as state variables, this function returns the parameters
of the stochastic version of that model

## Usage

``` r
parameterConversion(unit, reactants, products, kinetic.law)
```

## Arguments

- unit:

  of the kinetic parameters (a character vector), named.

- reactants:

  a list of the stoichiometric constants of each reaction, a list of
  integer vectors

- products:

  a list of the stoichiometric constants of each reaction, a list of
  integer vectors

- kinetic.law:

  a character matrix with two columns: ,1 for forward reaction rates,
  and ,2 for backward reaction rates.

## Value

a vector of parameter conversion factors

## Details

For this function, it is important that "2 A -\> B" is not written as
"A + A -\> B" (this may be fixed later).

## Examples

``` r
m <- model_from_tsv(uqsa_examples("AKAP79"))
#> Error in uqsa_examples("AKAP79"): could not find function "uqsa_examples"
r <- lapply(strsplit(m$Reaction$reactants,"+",fixed=TRUE),trimws)
#> Error in eval(expr, envir, enclos): object 'm' not found
p <- lapply(strsplit(m$Reaction$products ,"+",fixed=TRUE),trimws)
#> Error in eval(expr, envir, enclos): object 'm' not found
k <- ifelse(grepl("-",m$Reaction$kinetic.law,fixed=TRUE),m$Reaction$kinetic.law,paste0(m$Reaction$kinetic.law," - 0"))
#> Error in eval(expr, envir, enclos): object 'm' not found
f <- parameterConversion(u,stoichiometry(r),stoichiometry(p),k)
#> Error in eval(expr, envir, enclos): object 'u' not found
print(f)
#> Error in eval(expr, envir, enclos): object 'f' not found
```
