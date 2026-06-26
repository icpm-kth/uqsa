# Tabular Model File

``` r
library(uqsa)
library(errors)
```

The models we use in this package are written in a very loose tabular
form: the model is a collection of TSV files (tab-separated-values).
Tabs are the best separator in our opinion (as compared to commas or
semi-colon separated files), because it *can* be entered with a text
editor into the file, by hand, but is quite hard to enter as the content
of a field in a spreadsheet program. So, TABs usually don’t appear
coincidntally in the files, but can be typed on purpose. Furthermore:

- TSV files render nicely on sites such as GitHub
- TSV files are text files and thus work well with git
- TSV files do not require a library to read
- Column-names can contain spaces, which is very human friendly

We interpret the files partially using functions such as `grep` (in R),
`tolower`, with settings such as `ignore.case` where possible, so that
various different spellings are acceptable, e.g. “stdv”, “sd”, or
“standard deviation” when describing a distribution.

When reporting data values, we use parenthesised (concise) error
notation, e.g.: `1.23(1)` which is intrepreted as `1.23` with an
uncertainty of `0.01`, more examples:

|         string |     value | standard-error |
|---------------:|----------:|:---------------|
|      `1.23(1)` |      1.23 | 0.01           |
|    `1.23(1)e2` |       123 | 1              |
| `-1.234(1)e-6` | -1.234e-6 | 1e-9           |
|    `1.23;0.01` |      1.23 | 0.01           |
|    `1.23±0.01` |      1.23 | 0.01           |

The last two rows represent a fall-back notation with redundand zeros,
which you can use if you really don’t like parenthesised location, the
parser automatically tries to find some separator if there are no
parentheses in the string:

``` r
x <- parse_concise(
    c(
        "1.23(1)","1.23(1)e2","-1.234(1)e-6", # concise notation
        "1.23;0.01", "1.23±0.01"              # fall-back
    )
)
print(class(x))
#> [1] "errors"
print(as.data.frame(x)) # this looks better as a data-frame
#>              x
#> 1      1.23(1)
#> 2       123(1)
#> 3 -1.234(1)e-6
#> 4      1.23(1)
#> 5      1.23(1)
```

## Components of Models

All quantities in the model have a unit. If your model is initially
formulated using concentrations and reaction rates (kinetic laws), then
you can still simulate it using the Gillespie solver. We convert all
quantities with concentrations in the unit to particle counts
automatically. So, feel free to formulate the model using
concentrations.

The role of each file is determined from the file’s name:

- **Constant.tsv**, list of constants that are not subjetc to
  optimization, calibration, or any other investigation
- **Input.tsv**, known parameters that can change between different
  simulations
- **Parameter.tsv**, possibly unknown parameters (e.g. ranges), subject
  to fitting/sampling/etc.
- **Expression.tsv**
- **Compound.tsv**, list of reacting compounds, e.g. `Ca`
- **Reaction.tsv**, list of reactions
- **Output.tsv**, algebraic functions that express a measureable value
  (or very close to measureable),
  - e.g. the total amount of something: `A + AB + AC`, the total amount
    of bound and free `A`
- **Experiment.tsv**, list of experiments

The first column in each TSV table contains some kind of unique *id*,
the id of each *experiment* (row) in the experiment table
(Experiment.tsv) is also used as the name of the TSV file that contains
the data for this experiment.

## Data

The TSV files correspond pretty directly to R `data.frame` objects, the
import function `model_from_tsv` returns a list of data frames, named
like the file. The model’s name is taken from the directory (`dirname`)
the models are stored in. This name is stored in the `comment`
attribute.

There are two main cases we want to distinguish in the same way that
fungi and plants are distinct:

1.  Experimental data that corresponds fairly well to an output function
    - there is an output function that can be compared to the data,
      perhaps up to some scaling constant
    - the data column is labelled exactly like the output function is
      named, e.g.: `AMPA_OUTPUT`
2.  Experimental data that has a very complex relationship to the model
    - raw data, where several columns together can be used to make a
      comparison with something in the model

### One-to-One Correspondence

The first case, with a one-to-one correspondence between output
functions and data columns can be automatically parsed using the
`experiments` function. This function reads several of the tables and
determines the initial state, the input, the time-line for any given
experiment and creates a data-matrix (`data`), with standard-errors
using the `errors` package. Here, you can still create your own custom
likelihood function that performs some kind of complex normalization to
evaluate a simulation. It is not difficult to have a control experiment
that contains a value we need to normalize with. The user supplied
likelihood may perform such normalizations.

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
print(m$Reaction[,c(2,3,4)])           # an example table, the reactions
#>           reactants arw       products
#> r51           Rii_C <=>         RiiP_C
#> r14        RiiP + C <=>         RiiP_C
#> r12   RiiP_C + cAMP <=>    RiiP_C_cAMP
#> r43     cAMP + RiiP <=>      RiiP_cAMP
#> r23   RiiP_cAMP + C <=>    RiiP_C_cAMP
#> r78      cAMP + Rii <=>       Rii_cAMP
#> r56    Rii_C + cAMP <=>     Rii_C_cAMP
#> r76    Rii_cAMP + C <=>     Rii_C_cAMP
#> r62      Rii_C_cAMP <=>    RiiP_C_cAMP
#> r58         Rii + C <=>          Rii_C
#> r44      RiiP + CaN <=>       RiiP_CaN
#> r33 CaN + RiiP_cAMP <=>  RiiP_cAMP_CaN
#> r48        RiiP_CaN <=>      Rii + CaN
#> r37   RiiP_cAMP_CaN <=> CaN + Rii_cAMP
#> r1        C + AKAR4 <=>        AKAR4_C
#> r2          AKAR4_C <=>     AKAR4p + C
print(head(m$Experiment))   # the list of experiments
#>                          cAMP CaN b_AKAP  t0        type event
#> EX11____0nM__TRUE___TRUE  0.0 1.5      1 -30 time series      
#> EX12____0nM__TRUE__FALSE  0.0 1.5      0 -30 time series      
#> EX13____0nM_FALSE__FALSE  0.0 0.0      0 -30 time series      
#> EX21__100nM__TRUE___TRUE  0.1 1.5      1 -30 time series      
#> EX22__100nM__TRUE__FALSE  0.1 1.5      0 -30 time series      
#> EX23__100nM_FALSE__FALSE  0.1 0.0      0 -30 time series      
#>                                                                  comment
#> EX11____0nM__TRUE___TRUE Figure 3 in https://doi.org/10.7554/eLife.68164
#> EX12____0nM__TRUE__FALSE Figure 3 in https://doi.org/10.7554/eLife.68164
#> EX13____0nM_FALSE__FALSE Figure 3 in https://doi.org/10.7554/eLife.68164
#> EX21__100nM__TRUE___TRUE Figure 3 in https://doi.org/10.7554/eLife.68164
#> EX22__100nM__TRUE__FALSE Figure 3 in https://doi.org/10.7554/eLife.68164
#> EX23__100nM_FALSE__FALSE Figure 3 in https://doi.org/10.7554/eLife.68164
```

This is how the output function corresponds to the data-label,
**AKAR4pOUT**:

``` r
print(m$Output)
#>                        formula unit
#> AKAR4pOUT (AKAR4p*5)*71.67+100   µM
rn <- rownames(m$Experiment)
print(head(m[[rn[18]]]))    # experiment 18, data found by name
#>           time    AKAR4pOUT
#> E0301T001  -15  9.941(71)E1
#> E0301T002  -10   1.006(4)E2
#> E0301T003   -5 1.0194(92)E2
#> E0301T004    0         <NA>
#> E0301T005    5 1.0120(28)E2
#> E0301T006   10 1.0204(56)E2
```

So, the data-table (called like the row in the experiments table), has a
column that corresponds to the output function. This is how we link the
two together and know that the measured values have to somehow
correspond to the output function `AKAR4pOUT`.

``` r
x <- experiments(m)
print(x[[18]]$data[,seq(8),drop=FALSE])     # a sub-set of the data-matrix
#> Errors: 0.71 0.40 0.92   NA 0.28 ...
#>            [,1]  [,2]   [,3] [,4]  [,5]   [,6]   [,7]  [,8]
#> AKAR4pOUT 99.41 100.6 101.94   NA 101.2 102.04 102.35 103.8
```

The default likelihood function will take the squared difference between
the simulated values of `AKAR4pOUT` (the function) and the measured
values in that column.

### Very Indirect Data

In the second big case of very *cryptic* data, you should not rely on
the `data` matrix returned by the `experiments(m)` function and instead
use the `measurements` field, which is a data frame, exactly as it was
written in the TSV file, of just use the data data frame in `m`:

``` r
print(head(m[[rn[18]]],12))         # these two should be the same
#>           time    AKAR4pOUT
#> E0301T001  -15  9.941(71)E1
#> E0301T002  -10   1.006(4)E2
#> E0301T003   -5 1.0194(92)E2
#> E0301T004    0         <NA>
#> E0301T005    5 1.0120(28)E2
#> E0301T006   10 1.0204(56)E2
#> E0301T007   15 1.0235(43)E2
#> E0301T008   20   1.038(4)E2
#> E0301T009   25  1.055(11)E2
#> E0301T010   30    1.06(1)E2
#> E0301T011   35   1.106(5)E2
#> E0301T012   40  1.110(13)E2
print(head(x[[18]]$measurement,12)) # (a 12 row subset)
#>             time AKAR4pOUT
#> E0301T001 -15(0)   99.4(7)
#> E0301T002 -10(0)  100.6(4)
#> E0301T003  -5(0)  101.9(9)
#> E0301T004   0(0)    NA(NA)
#> E0301T005   5(0)  101.2(3)
#> E0301T006  10(0)  102.0(6)
#> E0301T007  15(0)  102.4(4)
#> E0301T008  20(0)  103.8(4)
#> E0301T009  25(0)    106(1)
#> E0301T010  30(0)    106(1)
#> E0301T011  35(0)  110.6(5)
#> E0301T012  40(0)    111(1)
```

The only difference is that all numbers were parsed by `parse_concise`
to resolve all parenthesised standard-errors (in `x`), in `m` the values
are raw (either stings if they contain parentheses, or numbers if R
could successfully coerce them.

So, `m` can be used directly, if you have a very complex data-case and
thus a very complex likelihood.
