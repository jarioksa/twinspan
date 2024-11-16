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

## Installation

You can install the current development version from GitHub using
**devtools** package:
```r
devtools::install_github("jarioksa/twinspan")
```
You need development tools to install a source package. In particular,
you need **C** and **Fortran** compilers.

If you cannot install a source package, you can install **twinspan** from R-Universe:

```r
install.packages('twinspan', repos = c('https://jarioksa.r-universe.dev', 'https://cloud.r-project.org'))
```

  
## What you can do with twinspan?

The basic command to run twinspan is – unsurprisingly – `twinspan`:
```r
> library(twinspan)
> library(natto) # for data 
> data(spurn)
> tw <- twinspan(spurn, cutlevels = 0:6)
```
TWINSPAN uses basically binary data, and quantitative data are broken into 
pseudospecies by species abundances using argument `cutlevels`. The default levels
are `c(0, 2, 5, 10, 20)` are suitable for percentage values, but the `spurn` data are cover classes.
Too keep the original class values, we list all cover classes as cut levels. Any species present
in the data will be marked as `species1` where the last digit shows the level of pseudospecies
level. With the present setting, `species3` means that the species occurs at cover class
value 3 or higher, and there will also be `species1` and `species2`. 

```r
> twindat <- twinsform(spurn, cutlevels=0:6)
> colnames(twindat)
  [1] "Elaerham1"  "Jacovulg1"  "Soladulc1"  "Rubufrut1"  "UrtidioiA1"
  [6] "Rumecris1"  "ClayperfB1" "StelmediB1" "FestrubrC1" "ElymrepeC1"
... cut ...
[111] "ElymrepeC4" "AmmoarenC4" "Epilangu4"  "Elaerham5"  "UrtidioiA5"
[116] "ClayperfB5" "StelmediB5" "FestrubrC5" "ElymrepeC5" "AmmoarenC5"
[121] "Elaerham6"  "ClayperfB6" "AmmoarenC6" "Elaerham7" 
```
The original data of 40 species are extended to a matrix of 124 pseudospecies.

Use `summary` to see the classification:
```r
> summary(tw)
1) eig=0.56:  +ElymrepeC1 < 1
  2) eig=0.399:  -UrtidioiA1 < 0
    4) eig=0.355:  -Soncarve1 < 0
      8) N=1: A19 
      9) N=4: A1 A2 A3 A20 
    5) eig=0.263:  +Inulcony1 < 1
      10) eig=0.229:  -Rumecris1 < 0
        20) N=1: B16 
        21) N=4: B4 B5 B10 B12 
      11) N=2: B17 B18 
  3) eig=0.23:  +AmmoarenC5 +Jacovulg1 -Galiveru1 < 1
    6) N=3: C6 C8 X14 
    7) eig=0.226:  -Calysold1 < 0
      14) N=3: C7 C9 C15 
      15) N=2: C11 C13 
```
`twinspan` is divisive: it splits data into two parts at each step, and these steps are
described here. The splits are based on the first correspondence analysis axis of the
current subset which is still further polished to make the dichotomy clearer. The first
split is made at eigenvalue 0.56, and the pseudospecies best indicating this division
are is `ElymrepeC1` (_Elymus repens_ at class value 1). It is an indicator of "positive"
(or right or second). The pseudospecies with minus sign indicate "negative" (or left or
first) group. The indicator species are summed up for every quadrat (or sampling unit: quadrat is
the term used in TWINSPAN). The last number after `<` gives the threshold score for
positive group: if the indicator score is less than 1, the quadrat is in the negative
group (2), and if it is 1 (or higher), the quadrat is in the positive group (3).
These groups are again divided with new correspondence analysis, and from group 2 you go 
either to 4 (negative) or 5 (positive). With default settings, groups smaller than 5
items or deeper than 7 levels of divisions are not divided. For these final groups,
`summary` gives the size (`N`) and lists the names of the member quadrats.

You can extract the classification of each quadrat with `cut`:
```r
> cut(tw)
 [1]  9  9  9 21 21  6 14  6 14 21 15 21 15  6 14 20 11 11  8  9
> cut(tw, level=2)
 [1] 4 4 4 5 5 6 7 6 7 5 7 5 7 6 7 5 5 5 4 4
```
You can also predict the membership of quadrats based on the indicator pseudospecies
and threshold score. This can also be done with argument `newdata` using data set that
contains same species, but is not used in developing the classification.
```r
> predict(tw, level=2)
 [1] 4 4 4 5 5 6 7 6 7 5 7 5 7 6 7 5 5 5 4 4
```
TWINSPAN classification is based
on the polished ordination, and the indicator pseudospecies only *indicate* this
division, and `predict` based on pseudospecies can give different classificatin than
the actual TWINSPAN that split the data by polished axis of correspondence analysis.

TWINSPAN stands for *two-way* indicator species analysis, and in addition to quadrat
classification it also classifies the species (not the pseudospecies):
```r
> summary(tw, "species")
1) eig=0.904
  2) eig=0.499
    4) eig=0.257
      8) eig=0.1
        16) eig=0.041
          32) N=3: ClayperfB StelmediB GeasfornB 
          33) N=4: CerafontB CirsvulgB Heraspho CardhirsB 
        17) N=2: Inulcony Bracruta 
      9) eig=0.217
        18) eig=0.109
          36) eig=0.024
            72) N=2: Hyporadi Arrhelat 
            73) N=4: UrtidioiA Soncaspe EurhpraeA Lophhete 
          37) N=1: Sambnigr 
        19) N=4: Soladulc Rubufrut Epilangu Hypncupr 
    5) N=3: Elaerham Jacovulg Rumecris 
  3) eig=0.212
    6) N=4: Soncarve Calysold Agrostol Verocham 
    7) eig=0.115
      14) eig=0.046
        28) eig=0.036
          56) N=2: FestrubrC BracalbiC 
          57) N=3: ElymrepeC AmmoarenC PoapratC 
        29) N=3: RanubulbC PlanlancC Cladranf 
      15) N=5: Ononspin Galiveru Bryuincl Syntrura Bovinigr 
```
Species classification is based on correspondence analysis where species are weighted by their
indicator potential for the quadrat classification. You can extract the classification vector
with `cut(tw, "species")`. 

The data can be tabulated with:
```r
> twintable(tw)
                                      
                  00000000000011111111
                  00000111111100011111
                  011110000011   00011
                       01111          
                                      
                  A   AB  BBBB  X  CCC
                  1AAA21BB1111CC1CC111
                  91230645027868479513
 00000  ClayperfB -----6566443--------
 00000  StelmediB -----2322543--------
 00000  GeasfornB ------22-222--------
 00001  CerafontB -----2--2233--------
 00001  CirsvulgB -----2----2---------
 00001   Heraspho ----------2---------
 00001  CardhirsB ----------22--------
 0001    Inulcony ----2-----23--------
 0001    Bracruta --23--2-2222--------
 001000  Hyporadi 2-------------------
 001000  Arrhelat 3-------------------
 001001 UrtidioiA 25334---------------
 001001  Soncaspe ----2---------------
 001001 EurhpraeA 33233---------------
 001001  Lophhete ----2---------------
 00101   Sambnigr ----22--------------
 0011    Soladulc 2442332232-22-22--2-
 0011    Rubufrut 223-4223222---2-----
 0011    Epilangu 4---22--23----2-----
 0011    Hypncupr -222---22-----2-----
 01      Elaerham 77777777777756644555
 01      Jacovulg 33232222---2--222322
 01      Rumecris 322222-----2222-2-22
 10      Soncarve 2-------2---222-223-
 10      Calysold --------22-2-33222--
 10      Agrostol ----2-----2--22-22--
 10      Verocham 2-------2-----2-2---
 11000  FestrubrC -------2----33322345
 11000  BracalbiC ------------2-2-3-2-
 11001  ElymrepeC ------------22445222
 11001  AmmoarenC ------------44262653
 11001   PoapratC -------------22-2232
 1101   RanubulbC --------------2-2-22
 1101   PlanlancC --------------2-2323
 1101    Cladranf -------------------2
 111     Ononspin ------------3-------
 111     Galiveru ------------2-2-----
 111     Bryuincl ------------2-2-----
 111     Syntrura -------------22-----
 111     Bovinigr --------------2-----
20 sites, 40 species
```
The strings of `0` and `1` in front of the species name and above quadrat name (or number)
give the steps of division. The numeric values in the table are the pseudospecies values
of the analysis.


### References

Hill, M.O. (1979) _TWINSPAN: A FORTRAN program for arranging multivariate
data in an ordered two-way table by classification of the individuals and
attributes_. Ecology and Systematics, Cornell University, Ithaca, NY.
