plot.clusboot <-
function (x, col=NULL, ...)
{
  Dmat <- 1-x$proportions

  stress.func <- function (y, delta)
    {
       Y <- matrix (y, ncol=2)
       dd <- dist(Y)
       sum((dd-delta)^2)/sum(dd^2)
    }
  Y <- cmdscale(Dmat)
  y <- optim(as.vector(Y), stress.func, delta=as.dist(Dmat))$par
  Y <- matrix(y, ncol=2)

  plot (Y[,1], Y[,2], asp=1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")

  cluster.vec <- x$clustering
  k <- nlevels(factor(cluster.vec))
  if (missing(col)) col <- rainbow(k)
  if (length(col)<k) col <- rep(col,k)

  for (i in 1:k)
    points(Y[cluster.vec==i,1], Y[cluster.vec==i,2], col=col[i], ...)

  Y
}
