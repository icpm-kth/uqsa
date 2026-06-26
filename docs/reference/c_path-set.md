# Add information about the model's C code

Adds the location of the model's C code (a file). The model is typically
a list of named numeric and named character vectors, which describe the
(interpreted) model.

## Usage

``` r
c_path(o) <- value
```

## Arguments

- o:

  the ODE , or CME model

- value:

  the path to the compiled model

## Value

modified o, with information about compiled code m \<-
model_from_tsv(uqsa_example("AKAR4")) o \<- as_ode(m) c_path(o) \<-
write_c_code(generate_code(o)) so_path(o) \<- shlib(o) print(o)
