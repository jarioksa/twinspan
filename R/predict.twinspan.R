### predict classification using officia indicator pseudospecies

`predict.twinspan` <-
    function(object, newdata, ...)
{
    if (missing(newdata))
        stop("needs data to predict...")
    inds <- object$quadrat$indicators
    poslim <- object$quadrat$positivelimit
    indlab <- object$quadrat$indlabels
    cuts <- object$cutlevels
    ## handle newdata
    newdata <- twinsform(newdata, cuts, downweight = FALSE)
    ## one quadrat is turned into 1-column matrix: transpose
    if (NCOL(newdata) == 1)
        newdata <- t(newdata)
    m <- colnames(newdata) %in% indlab
    newdata <- newdata[,m, drop=FALSE]
    pred <- numeric(nrow(newdata))
    prow <- numeric(length(indlab))
    names(prow) <- indlab
    ## cycles over rows of newdata
    for(i in seq_len(nrow(newdata))) {
        prow[] <- 0
        prow[colnames(newdata)] <- newdata[i,]
        k <- 1
        repeat{
            ind <- inds[,k]
            if (all(ind == 0)) break
            ind <- ind[ind != 0]
            score <- sum(sign(ind) * prow[abs(ind)])
            if (score < poslim[k]) {
                k <- 2*k
            } else {
                k <- 2*k + 1
            }
            pred[i] <- k
            if (k > ncol(inds)) break
        }
    }
    pred
}
