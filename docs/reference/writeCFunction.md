# Write C Function

This function writes a C function according to a typical template.

## Usage

``` r
writeCFunction(
  prefix,
  fName,
  defArgs = c("double t", "const double y_[]"),
  retValue = NULL,
  defs = NULL,
  values = NULL,
  body = NULL,
  otherArgs = "void *par",
  init0 = TRUE
)
```

## Arguments

- prefix:

  function's prefix (name of model)

- retValue:

  name of return value (attr: memset?)

- defs:

  optional expressions, parameters, etc.

- values:

  a character vecotor with values or mathematic formulae to return.

- body:

  a character vector of additional content for the body of the function

- otherArgs:

  additional arguments

- fNname:

  function's name
