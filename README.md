# twinspan
**R** package for Two-Way Indicator Species Analysis (Hill 1979).

Two-Way Indicator Species Analysis was developed to classify
community data tables. It was supposed to work in the same 
way as a traditional community ecologist in arranging a
community data table.

TWINSPAN is available as a self-standing computer program that
can be compiled to work on many platforms. This R package uses
the same Fortran code, but allows using TWINSPAN from **R** 
together with other **R** functions for community ecology and
statistics.

The design philosophy of TWINSPAN is completely different from
well-behaved **R** functions. TWINSPAN is a traditional console
program that runs through a process, and prints its results as
it advances. In **R**, the function should work silently, and
return the result for further analysis. The code needs a thorough
re-design to be used in **R**.

## Version History

- **0.1:** Development version unsuitable for any work.
- **0.2:** Marginally usable: performs SU (quadrat) classification and has 
  support function `cut` to find SU class membership vector for any level
  of hierarchy. No species classification, no diagnostics, no information
  on indicator species.
- **0.3:** Adds species classification. 
- **0.4:** Returns eigenvalues and indicator species for each division.
  However, most support functions are still missing, and the result object
  must be inspected manually.
- **0.5:** Basically done and usable: output contains all basic
  objects that are needed. However, most support functions are
  missing, and these elements must be accessed and handled
  manually. The next section describes the structure of the result
  object and its handling.
- **0.6:** Has now most support functions: `summary` gives the division
  history and lists the signed indicator species, `cut` returns quadrat
  or species classification at any level, `twintable` prints the 
  classified community table, `as.dendrogram` displays the divisions as
  a denrogram that can be plotted and handled with standdard dendrogram 
  tools, `predict` uses indicator species to allocate quadrats to established
  classes also with `newdata`. In addition, `twinsform` transform data so that
  its correspondence analysis is similar to the one in `twinspan`. However,
  the function is not tested with all non-default options and documentation is
  still rudimentary.
  
## Structure of the Result Object

The structure of the result object is:
```r
List of 3
 $ call   : language twinspan(x = varespec)
 $ quadrat:List of 7
  ..$ iclass       : int [1:24] 39 31 14 6 30 6 31 31 6 5 ...
  ..$ eig          : num [1:63] 0.179 0.147 0.203 0.163 0 ...
  ..$ indicators   : int [1:7, 1:63] -105 91 0 0 0 0 0 45 76 0 ...
  ..$ positivelimit: int [1:63] 1 2 1 1 0 0 0 0 0 0 ...
  ..$ labels       : chr [1:24] "18" "15" "24" "27" ...
  ..$ indlabels    : chr [1:107] "Callvulg1" "Empenigr1" "Rhodtome1" "Vaccmyrt1" ...
  ..$ index        : int [1:24] 19 20 21 22 18 14 15 1 16 17 ...
 $ species:List of 4
  ..$ iclass: int [1:44] 10 64 30 30 64 64 31 14 33 11 ...
  ..$ eig   : num [1:63] 0.543 0.41 0.405 0.122 0.272 ...
  ..$ labels: chr [1:44] "Callvulg" "Empenigr" "Rhodtome" "Vaccmyrt" ...
  ..$ index : int [1:44] 2 5 6 17 19 27 28 29 30 43 ...
 - attr(*, "class")= chr "twinspan"
```

The example is from version 0.5, and some changes can be expected.
However, basically the structure is similar as planned in advance. The
main results are given in elements `quadrat` and `species` which give
the basic information of the two-way classification. These have similar
elements:

- `iclass`: The final classification vector. This gives the ID numbers
  of the deepest level of classification, but function `cut` returns
  the classification vector at any level of division. Function 
  `as.dendrogram` will turn this into a standard **R** dendrogam
  object that can be plotted and further manipulated.
- `eig`: Eigenvalue of the division. This is 0 for divisions that were
  skipped because group size was too small. These can be accessed with
  function `eigenvals`.
- `indicators`: Signed indices of indicator pseudospecies. The
  absolute value gives the index, and the sign tells if the
  pseudospecies is an indicator of the negative or positive
  group. This is a matrix where each column lists the indicators of
  the corresponding division, and if division was skipped, the column
  is zero.
- `positivelimit`: The lowest indicator value for the positive group;
  if the indicator score is less, the item goes to the negative group
  ("negative" and "positive" are just conventional names used in
  `twinspan` to separate the branches of dichotomy).
- `labels`: Names of the elements (quadrats or species).
- `indlabels`: Labels for the pseudospecies. These are made by
  appending the cutlevel number to species name. The absolute values
  of `indicators` index these labels. In this case, the first division
  has two indices `-105` and `91`: `indlabels[105]` is `Cladrang5`
  (which is a negative indicator at cutlevel 5), and `indlabels[91]`
  is `Pleuschr` (which is a positive indicator at cutlevel 4 or
  higher).
- `index`: Index that will order rows and columns according to the
  classification as returned from the Fortran code. The data can be
  tabulated by using these as arguments `site.ind` and `sp.ind` in
  **vegan** functions `vegemite` and `tabasco`.
  
### References

Hill, M.O. (1979) _TWINSPAN: A FORTRAN program for arranging multivariate
data in an ordered two-way table by classification of the individuals and
attributes_. Ecology and Systematics, Cornell University, Ithaca, NY.
