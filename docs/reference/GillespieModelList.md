# Create a list of reactions for Stochastic Simulations

This function takes a series of SBtab tables, as returned by
[`SBtabVFGEN::sbtab_from_tsv()`](https://rdrr.io/pkg/SBtabVFGEN/man/sbtab_from_tsv.html)
and creates GillespieSSA2 reactions from them. Reactions arfe made
pairwise, as forward and backward reaction pairs. If a backward reaction
doesn't exist, the list item is NULL. A valid set of reactions can be
obtained with `!is.null(reactions)`

## Usage

``` r
GillespieModelList(SBtab, LV = NULL)
```

## Arguments

- SBtab:

  a series of tables as returned by
  [`sbtab_from_tsv()`](https://rdrr.io/pkg/SBtabVFGEN/man/sbtab_from_tsv.html)

- LV:

  is the product of Avogadro's constant L and the system's volume V in
  litres; if unspecified this information is retrieved from the SBtab
  files, if missing we assume 1µm³ of volume (the approximate sizes of
  bacteria or synapses)

## Value

a list of reactions

## Examples

``` r
 # model.tsv <- dir(pattern="[.]tsv$")
 # model.sbtab <- SBtabVFGEN::sbtab_from_tsv(model.tsv)
 # reactions <- makeGillespieModel(model.sbtab)
 # l <- is.null(reactions)
 # model.ssa2 <- reactions[!l]
```
