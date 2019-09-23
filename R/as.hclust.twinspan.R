### Quadrat divisions without terminal items (quadrats) as an hclust
### tree. Hclust trees can not have polytomies, and quadrats must be
### omitted. There already is as.dendrogram, but try to see if this
### can be developed into a decision tree.

`as.hclust.twinspan` <-
    function(x, ...)
{
    class <- cut(x)
    n <- max(class)
    state <- character(n)
    state[which(x$quadrat$eig > 0)] <- "div"
    state[unique(class)] <- "clust"
    nclus <- sum(state == "clust")
    line <- numeric(n)
    merge <- matrix(0, nclus-1, 2)
    height <- numeric(nclus-1)
    pow2 <- 2^(1:15)
    maxh <- sum(max(class) >= pow2) + 1
    nro <- nclus
    now <- 0
    for(i in seq(n-1, 2, by=-2)) {
        if (state[i] == "")
            next
        now <- now + 1
        height[now] <- maxh - sum(i >= pow2)
        for(j in 2:1) {
            line[i+j-1] <- now
            if(state[i+j-1] == "clust") {
                merge[now,j] <- -nro
                nro <-  nro - 1
            } else {
                merge[now,j] <- line[(i + j - 1) * 2]
            }
        }
    }
    labels <- which(state=="clust")
    ind <- x$quadrat$index
    order <- order(tapply(order(ind), class, mean))
    out <- list(merge = merge, labels = labels, height = height, order = order,
                method = "twinspan")
    class(out) <- "hclust"
    out
}
