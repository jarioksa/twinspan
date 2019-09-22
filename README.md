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
  a denrogram that can be plotted and handled with standard dendrogram 
  tools, `predict` uses indicator species to allocate quadrats to established
  classes also with `newdata`. In addition, `twinsform` transforms data so that
  its correspondence analysis is similar to the one in `twinspan`. However,
  the function is not tested with all non-default options and documentation is
  still rudimentary.
  
## What you can do with twinspan?

The basic command to run twinspan is – unsurprisingly – `twinspan`:
```r
> library(twinspan)
> library(vegan) # for data 
> data(varespec)
> tw <- twinspan(varespec)
```
This will run the analysis with default options, including the default cut levels.
TWINSPAN uses basically binary data, and quantitative data are broken into 
pseudospecies by species abundances using argument `cutlevels`. The default levels
are `c(0, 2, 5, 10, 20)`. Any species present in the data will be marked as `species1`
where the last digit shows the level of pseudospecies level. If that species is present
at abundance 7, TWINSPAN will also generate `species2` and  `species3`, because these
limits are exceeded. To see how data are transformed, you use command `twinsform`
that transforms the data similarly as `twinspan`:
```r
> twindat <- twinsform(varespec)
> colnames(twindat)
  [1] "Callvulg1" "Empenigr1" "Rhodtome1" "Vaccmyrt1" "Vaccviti1" "Pinusylv1"
  [7] "Descflex1" "Betupube1" "Vacculig1" "Diphcomp1" "Dicrsp1"   "Dicrfusc1"
 [13] "Dicrpoly1" "Hylosple1" "Pleuschr1" "Polypili1" "Polyjuni1" "Polycomm1"
 [19] "Pohlnuta1" "Ptilcili1" "Barbhatc1" "Cladarbu1" "Cladrang1" "Cladstel1"
 [25] "Cladunci1" "Cladcocc1" "Cladcorn1" "Cladgrac1" "Cladfimb1" "Cladcris1"
 [31] "Cladchlo1" "Cladbotr1" "Cladamau1" "Cladsp1"   "Cetreric1" "Cetrisla1"
 [37] "Flavniva1" "Nepharct1" "Stersp1"   "Peltapht1" "Icmaeric1" "Cladcerv1"
 [43] "Claddefo1" "Cladphyl1" "Callvulg2" "Empenigr2" "Rhodtome2" "Vaccmyrt2"
 [49] "Vaccviti2" "Descflex2" "Vacculig2" "Diphcomp2" "Dicrsp2"   "Dicrfusc2"
 [55] "Dicrpoly2" "Hylosple2" "Pleuschr2" "Polyjuni2" "Ptilcili2" "Barbhatc2"
 [61] "Cladarbu2" "Cladrang2" "Cladstel2" "Cladunci2" "Flavniva2" "Nepharct2"
 [67] "Stersp2"   "Callvulg3" "Empenigr3" "Vaccmyrt3" "Vaccviti3" "Vacculig3"
 [73] "Dicrsp3"   "Dicrfusc3" "Hylosple3" "Pleuschr3" "Polyjuni3" "Ptilcili3"
 [79] "Cladarbu3" "Cladrang3" "Cladstel3" "Cladunci3" "Flavniva3" "Stersp3"  
 [85] "Callvulg4" "Empenigr4" "Vaccmyrt4" "Vaccviti4" "Dicrsp4"   "Dicrfusc4"
 [91] "Pleuschr4" "Ptilcili4" "Cladarbu4" "Cladrang4" "Cladstel4" "Cladunci4"
 [97] "Flavniva4" "Stersp4"   "Callvulg5" "Vaccviti5" "Dicrsp5"   "Dicrfusc5"
[103] "Pleuschr5" "Cladarbu5" "Cladrang5" "Cladstel5" "Cladunci5"
```
The original data of 44 species are extended to a matrix of 107 pseudospecies.
People often use the default cut levels which are rather good for data originally
expressed in cover percentages. If you have cover class data and you want to
preserve the original accuracy, you can give your data values. For instance,
for the `dune` data in **vegan** you may use `twinspan(dune, cutlevels=1:9)`.

Use `summary` to see the classification:
```r
> summary(tw)
1) eig=0.179:  -Cladrang5 +Pleuschr4 < 1
  2) eig=0.147:  +Callvulg2 +Pleuschr3 < 2
    4) eig=0.163:  +Cladarbu3 < 1
      8) N=4: 2 9 12 10 
      9) eig=0.182:  -Cetrisla1 < 0
        18) N=1: 4 
        19) eig=0.169:  +Callvulg1 < 1
          38) N=2: 7 5 
          39) N=3: 18 6 3 
    5) N=3: 13 14 11 
  3) eig=0.203:  +Cetreric1 +Cladarbu2 -Cladstel2 < 1
    6) N=4: 27 19 28 21 
    7) eig=0.161:  -Dicrsp2 < 0
      14) N=2: 24 25 
      15) eig=0.206:  +Callvulg1 < 1
        30) N=1: 23 
        31) N=4: 15 22 16 20
```
`twinspan` is divisive: it splits data into two parts at each step, and these steps are
described here. The splits are based on the first correspondence analysis axis of the
current subset which is still further polished to make the dichotomy clearer. The first
split is made at eigenvalue 0.179, and the pseudospecies best indicating this division
are `Cladrang5` and `Pleuschr4`. `Cladrang5` is an indicator of "negative" (or left or
first) group and `Pleuschr4` an indicator of "positive" (or right or second) group. 
These indicator species are summed up for every quadrat (or sampling unit: quadrat is
the term used in TWINSPAN). The last number after `<` gives the threshold score for
positive group: if the indicator score is less than 1, the quadrat is in the negative
group (2), and if it is 1 (or higher), the quadrat is in the positive group (3).
These groups are again divided with new correspondence analysis, and from group 2 you go 
either to 4 (negative) or 5 (positive). With default settings, groups smaller than 5
items or deeper than 7 levels of divisions are not divided.

You can extract the classification of each quadrat with `cut`:
```
> cut(tw)
 [1] 39 31 14  6 30  6 31 31  6  5  5 31 14 38 38 39 39 18  8  8  8  8  5  6
> cut(tw, level=2) # use classification at second level
 [1] 4 7 7 6 7 6 7 7 6 5 5 7 7 4 4 4 4 4 4 4 4 4 5 6
```
You can also predict the membership of quadrats based on the indicator pseudospecies
and threshold score. This can also be done with argument `newdata` using data set that
contains same species, but is not used in developing the classification.
```r
> predict(tw, level=2)
 [1] 4 7 7 6 7 6 7 7 6 5 5 7 7 4 4 4 4 4 4 4 4 4 5 4
```
Please note that the last quadrat was classified to second-level class 6, but it is
predicted to be in class 4. This can happen because TWINSPAN classification is based
on the polished ordination, and the indicator pseudospecies only *indicate* this
division. The analysis tries to make the concordance as good as possible, but it
does not always succeed. Traditionally these quadrats are called misclassifications.
In this case, the last quadrat was misclassified on the first step: the ordination 
put it into group 3, but indicator pseudospecies to group 2.

TWINSPAN stands for *two-way* indicator species analysis, and in addition to quadrat
classification it also classifies the species (not the pseudospecies):
```r
> summary(tw, "species")
1) eig=0.543
  2) eig=0.41
    4) eig=0.122
      8) eig=0.068
        16) eig=0.049
          32) eig=0.048
            64) N=10: Empenigr Vaccviti Pinusylv Polyjuni Pohlnuta Cladcorn Cladgrac Cladfimb Cladcris Claddefo 
            65) N=3: Cladunci Cetrisla Peltapht 
          33) N=3: Vacculig Polypili Cladsp 
        17) N=4: Cladarbu Cladcocc Cetreric Cladcerv 
      9) N=4: Cladrang Cladchlo Nepharct Stersp 
    5) eig=0.272
      10) N=4: Callvulg Cladstel Icmaeric Cladphyl 
      11) N=3: Diphcomp Cladamau Flavniva 
  3) eig=0.405
    6) N=3: Dicrfusc Dicrpoly Ptilcili 
    7) eig=0.231
      14) N=3: Betupube Dicrsp Cladbotr 
      15) eig=0.254
        30) N=3: Rhodtome Vaccmyrt Hylosple 
        31) N=4: Descflex Pleuschr Polycomm Barbhatc 
```
Species classification is based on correspondence analysis where species are weighted by their
indicator potential for the quadrat classification. You can extract the classification vector
with `cut(tw, "species")`. 

The data can be tabulated with:
```r
> twintable(tw)
                                         
                 000000000000011111111111
                 000000000011100001111111
                 0000111111       0011111
                     011111         01111
                      00111              
                                         
                   11   1  11121222221212
                 292047586334179814535260
 000000 Empenigr 324413143213143141341333
 000000 Vaccviti 354413244233444354454424
 000000 Pinusylv 11111-111-111-1111111111
 000000 Polyjuni 11-111111-111-1111321-11
 000000 Pohlnuta -11111-1-111111111-11111
 000000 Cladcorn 111111111111111111111111
 000000 Cladgrac 111-11111111111111111111
 000000 Cladfimb 111111111111-1111-111111
 000000 Cladcris 111111111111111111111111
 000000 Claddefo 111111111111-11111111111
 000001 Cladunci 111111112115211113112212
 000001 Cetrisla -1111-----111--1111-1-1-
 000001 Peltapht -------1--1-11----11----
 00001  Vacculig ----13-1-1---1-1--1--2-1
 00001  Polypili ------111---1-1----11--1
 00001    Cladsp -1--1--11-1--1----1--111
 0001   Cladarbu 112145455354313112334333
 0001   Cladcocc -1-1111111111-11--111111
 0001   Cetreric -1-111111111-----11-1111
 0001   Cladcerv 1---1-------------1-----
 001    Cladrang 524355555551532123233244
 001    Cladchlo 11-1-1--111---1-11-1-1--
 001    Nepharct -----1-1-1--------2----1
 001      Stersp ---1114111111-11-1111-1-
 010    Callvulg 1-112--111522----1--1221
 010    Cladstel 555551124541525-41-11111
 010    Icmaeric -----11----1----------1-
 010    Cladphyl -1-1-------11-----------
 011    Diphcomp -1---1-2-11-------1-----
 011    Cladamau ------11-1--------------
 011    Flavniva -1-1411111----1----1----
 10     Dicrfusc 11111111112412111-424551
 10     Dicrpoly -1-1-11-----1-1121-11-11
 10     Ptilcili 11-----1111-11114111111-
 110    Betupube ----------------1-1--1--
 110      Dicrsp --1-------1----1154-1111
 110    Cladbotr ----------1--11111-1---1
 1110   Rhodtome ----------1--2-12----1--
 1110   Vaccmyrt -1--------11-3244---132-
 1110   Hylosple -------------3-3--1----1
 1111   Descflex ----11-------2111-1--1--
 1111   Pleuschr 123111121133455525545544
 1111   Polycomm -----------1-11-1-1-----
 1111   Barbhatc ----------1--11-2--1----
  sites species 
     24      44 
```
The strings of `0` and `1` in front of the species name and above quadrat name (or number)
give the steps of division. The numeric values in the table are the pseudospecies values
of the analysis.

The `twinspan` classification can be extracted as a standard **R** `dendrogram`. The final
units contain many branches (species, quadrats), and it is best to use fan-like trees:
```r
plot(as.dendrogram(tw, "species"), type = "triangle")


### References

Hill, M.O. (1979) _TWINSPAN: A FORTRAN program for arranging multivariate
data in an ordered two-way table by classification of the individuals and
attributes_. Ecology and Systematics, Cornell University, Ithaca, NY.
