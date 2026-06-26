# scrnn returns a closure around gsl_odeiv2_CRNN()

the returned value is a function of a variable p that encodes the CRNN
in some way. Three user supplied functions are used to extract the three
components of a CRNN:

## Usage

``` r
scrnn(
  experiments,
  modelName,
  parMap = function(p) p$l,
  stoichiometry = function(p) p$nu,
  modifiers = function(p) p$m,
  method = 0,
  time.out = 1
)
```

## Arguments

- experiments:

  list of experiments (inputs are ignored).

- modelName:

  scalar string, can indicate a shared library with an attached comment
  attribute.

- parMap:

  (function) extracts kinetic rate coefficients from its argument.

- stoichiometry:

  (function) extracts the stoichiometry matrix from its argument.

- modifiers:

  (function) extracts the modifier matrix from its argument.

- method:

  (integer) integration method key (0:10) corresponds to these GSL
  methods: msbdf, msadams, bsimp, rk4imp, rk2imp, rk1imp, rk8pd, rkck,
  rkf45, rk4, rk2

- time.out:

  time limit for solution in seconds

## Value

closure that maps one argument (p) to simulation results (y).

## Details

- kinetic rate coefficients (in log-space): `l <- parMap(p)`

- stoichiometric matrix: `nu <- stoichiometry(p)`

- modifier matrix: `m <- modifiers(p)`

these three components (one numeric vector, and two matrices) are passed
to the simulation procedure. The vector l can be a matrix with M
columns. In that case, one simulation per column is performed. The
stoichiometry and modifiers remain unchanged throughout.

## Examples

``` r
# \donttest{
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  ex <- experiments(m)
  nu <- stoichiometric_matrix(m)
  l <- matrix(
    c(log(values(m$Parameter)),-1e6),
    2,2,
    byrow=TRUE,
    dimnames=list(rownames(m$Reaction),c("fwd","bwd"))
  )
  C <- CRNN(
    NCOL(nu),
    initialValues=values(m$Compound),
    funcValues=formulae(m$Output)
  )
  c.file <- tempfile("AKAR4",fileext=".c")
  cat(C,file=c.file,sep='\n')
  modelName <- "CRNN"
  comment(modelName) <- shlib(c.file)
  s <- scrnn(ex, modelName)
  p <- list(l=l,nu=nu,m=nu*0)
  y <- s(p)
  plot(ex,y)



#> NULL
# }
```
