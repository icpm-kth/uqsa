# uqsa 0.6.3

The way that we generate C code and how shared libraries are made has
changed. In the past we were aiming for a workflow, where a user would
have their model in a repository. This repository would have the
generated C code for their model as well, so that there was no need to
keep regenerating it. The same repository could have a copy of the
model in various formats such as `mod`, `sbml`, `vf`, etc.

However, these derived formats make it quite difficult to keep
everything synchronized and up to date. Therefore we consider it the
more common case that we frequently change the model TSV files and
then update the C code.

The model, once converted to an ODE or CME (for chemical master
equation mtype odels) now has a field for the C file storage location
and the `shlib` function uses this field to compile the model to a
shared library (same directory). The path to the `so` file is also
stored in the ODE or CME object.

Many other changes have been made to improve convenience, as well es
to fix bugs.
