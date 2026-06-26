# Retrieve information about compiled code

Returns the path of the shared library (.so file). The model is
typically a list of named arrays and matrices.

## Usage

``` r
so_path(o)
```

## Arguments

- o:

  the ODE, or CME model

## Value

modified o, with information about compiled code m \<-
model_from_tsv(uqsa_example("AKAR4")) o \<- as_ode(m) c_path(o) \<-
write_c_code(generate_code(o)) so_path(o) \<- shlib(o) print(so_path(o))
