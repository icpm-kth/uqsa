# model_from_tsv loads the content from a series of tsv files

The argument can be either a series of tsv file-names, or a directory
with tsv files. If it is a directory, all tsv files therein wil be used.

## Usage

``` r
model_from_tsv(src = ".")
```

## Arguments

- src:

  either a vector of files, or a directory with tsv files

## Value

a list of data.frames, one per file, named like the files.
