### TWINSPAN table

`twintable` <-
    function(object)
{
    i <- object$quadrat$index
    j <- object$species$index
    mat <- twin2mat(object)
    vegan::vegemite(mat[i,j], zero="-")
}
