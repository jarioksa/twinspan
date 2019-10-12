### sum of all eigenvalues in Correspondence Analysis (total Chi-squared)

`totalchi` <-
    function(x)
{
    x <- x/sum(x)
    rc <- outer(rowSums(x), colSums(x))
    x <- (x - rc)/sqrt(rc)
    sum(x^2)
}

### Evaluates the sum of all eigenvalues of Correspondence Analysis
### for all divisions and terminal groups > 1 member. Twinspan only
### evaluated the first eigenvalue of divisions

`twintotalchi` <-
    function(x, what = c("quadrat", "species"))
{
    what <- match.arg(what)
    chi <- numeric(2^(x$levelmax+1) - 1)
    for(lev in 0:x$levelmax) {
        ids <- cut(x, level = lev, what = what)
        tab <- table(ids)
        for (k in unique(ids))
            if(sum(ids==k) > 1) {
                z <- switch(what,
                    "quadrat" =
                        twin2stack(x, subset = ids==k, downweight = TRUE),
                    "species" =
                        twin2specstack(x, subset = ids==k, downweight = TRUE))
                chi[k] <- totalchi(z)
            }
    }
    chi
}

