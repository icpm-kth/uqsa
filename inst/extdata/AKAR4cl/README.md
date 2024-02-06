# AKAR4 cl

This is a version of the model that was built with conservation law
analysis turned on when processing the SBtab content to produce the
vector field file (`vf`). This has a downstream effect on the R and C
code.

The state variables `AKAR4` and `AKAR4_C` have been replaced by the equivalent expressions:

```R
AKAR4_C <- AKAR4_C_ConservedConst - (C)
AKAR4 <- AKAR4_ConservedConst - (AKAR4p-C)
```

where the constants are derived from the initial values specified in SBtab:

```txt
!!SBtab Document='AKAR4cl' TableName='Compound' SBtabVersion='1.0'
!ID      !InitialValue
AKAR4    0.2
AKAR4_C  0
AKAR4p   0
C        0
```

These constants are declared as *inputs* in the subsequent code. The
[`SBtabVFGEN::sbtab.data`](https://github.com/icpm-kth/SBtabVFGEN)
function will respect the conservation law analysis results and
redeclare initial *experiment-specific initial values* as
*experiment-specific inputs*. For this to work, the optional `conLaws`
argument must be provided like in the [runAKAR4cl](runAKAR4cl)
script. This model doesn't declare any other inputs, so afterwards the
input entries `experiments[[i]]$input` should contain only these
components:

```R
> experiments[[1]]$input
AKAR4_C_ConservedConst   AKAR4_ConservedConst
                   0.4                   -0.2
```

and all of them reduced to one matrix:

```R
> inputs<-Reduce(function(a,b) {rbind(a,b$input)},experiments,init=NULL)
     AKAR4_C_ConservedConst AKAR4_ConservedConst
[1,]                  0.400               -0.200
[2,]                  0.100                0.100
[3,]                  0.025                0.175
```

