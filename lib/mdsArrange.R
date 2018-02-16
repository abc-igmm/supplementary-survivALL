
mdsArrange <- function(d, isAlreadyScaled = FALSE){
    if(isAlreadyScaled == FALSE){
        dst <- dist(t(scale(d)))
        dN <- dimnames(d)[[2]]
        dst.m <- as.matrix(dst)
        dimnames(dst.m) <- list(dN, dN)
        dst.cmd <- cmdscale(dst.m, k=2)
    } else {
        dst <- dist(t(d))
        dN <- dimnames(d)[[2]]
        dst.m <- as.matrix(dst)
        dimnames(dst.m) <- list(dN, dN)
        dst.cmd <- cmdscale(dst.m, k=2)
    }

  data.frame(ids = row.names(dst.cmd), x = dst.cmd[,1], y = dst.cmd[,2])

}
