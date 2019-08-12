# twinspan
R package for Two-Way Indicator Species Analysis (Hill 1979)

Two-Way Indicator Species Analysis was developed to classify
community data tables. It was supposed to work in the same 
way as a traditional community ecologist in arranging a
community data table.

TWINSPAN is available as a self-standing computer program that
can be compiled to work on many platforms. This R package uses
the same Fortran code, but allows using TWINSPAN from **R** 
together with other **R** functions for community ecology and
statistics.

## Status of the Package

The design philosophy of TWINSPAN is completely different from
well-behaved **R** functions. TWINSPAN is a traditional console
program that runs through a process, and prints its results as
it advances. In **R**, the function should work silently, and
return the result for further analysis. The code needs a thorough
re-design to be used in **R**. Currently, the package can only
return the final classification for quadrats. In the future, similar
species classification will be available. The following features
need more work:

- Information on indicator species for each division.
- Diagnostics for indicator divisions.
- Eigenvalues associated with each division.
- Handling of species and SU names.
- Reporting misclassifications (_i.e._, indicator species give different
  results than actual classification).
- Support functions (plotting, printing, extracting classifications)
