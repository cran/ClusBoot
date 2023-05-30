# Sugnet Lubbe
# MuViSU (Centre for Multi-Dimensional Data Visualisation)
# Department of Statistics and Actuarial Science
# Stellenbosch University, June 2022
# ============================================================

# ---------------------------------------------------------------------------------
# clusboot
# ---------------------------------------------------------------------------------

#' Performs bootstrap on a cluster analysis output
#'
#' @param datmat a data matrix or distance object which will be the input to the clustering function
#' @param B number of bootstrap replicates
#' @param clustering.func the function which will perform the clustering and output a vector of cluster memberships
#' @param ... more arguments to be passed to the clustering function
#'
#' @return an object of type clusboot
#' @export
#'
#' @examples
#' clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)

clusboot <- function(datmat, B=1000, clustering.func=complete.linkage, ...)
{
  obj <- datmat
  if (data.class(datmat) != "matrix") datmat <- as.matrix(datmat)
  n <- nrow(datmat)
  if (is.null(rownames(datmat))) sample.names <- 1:n
  else sample.names <- rownames(datmat)
  boot.samples <- matrix(sample (1:n, size=n*B, replace=T), ncol=B)

  boot.out <- apply(boot.samples,2,function(x)
  { if (data.class(obj) == "dist") out <- clustering.func(stats::as.dist(datmat[x,x]), ...)
  else out <- clustering.func(datmat[x,], ...)
  clusD <- totD <- matrix (0, nrow=n, ncol=n, dimnames=list(sample.names,sample.names))

  boot.sample <- table(x)
  boot.names <- as.numeric(names(boot.sample))
  Dmat <- matrix (boot.sample,ncol=1) %*% matrix (boot.sample, nrow=1)
  totD[boot.names,boot.names] <- Dmat

  kk <- nlevels(factor(out))
  for (i in 1:kk)
  {
    clus.count <- table(x[out==i])
    clus.names <- as.numeric(names(clus.count))
    Dmat <- matrix (clus.count,ncol=1) %*% matrix (clus.count,nrow=1)
    clusD[clus.names,clus.names] <- Dmat
  }

  cbind(clusD, totD)
  })
  # boot.out is a matrix with B columns
  # the first nrow/2 is the n*n elements of clusD
  # the remaining nrow/2 is the n*n elements of totD

  boot.out <- apply (boot.out, 1, sum)
  clusD <- matrix (boot.out[1:(n*n)], nrow=n, ncol=n)
  totD <- matrix (boot.out[-(1:(n*n))], nrow=n, ncol=n)

  clusD <- clusD/totD
  dimnames(clusD) <- list(sample.names,sample.names)
  out <- clustering.func(obj, ...)
  output <- list(proportions=clusD[order(out),order(out)], clustering=out[order(out)])
  class(output) <- "clusboot"
  output
}

# ---------------------------------------------------------------------------------
# plot.clusboot
# ---------------------------------------------------------------------------------

#' MDS plot of similarities given by the proportion of bootstrap replicates where objects cluster together
#'
#' @param x an object of class clusboot
#' @param col single colour or a vector specifying a colour for each object
#' @param ... more arguments to be passed to `plot()`
#'
#' @return matrix of similarities (proportions)
#' @export
#'
#' @examples
#' out <- clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' plot(out)

plot.clusboot <- function(x, col, ...)
{
  Dmat <- 1-x$proportions

  stress.func <- function (y, delta)
  {
    Y <- matrix (y, ncol=2)
    dd <- stats::dist(Y)
    sum((dd-delta)^2)/sum(dd^2)
  }
  Y <- stats::cmdscale(Dmat)
  y <- stats::optim(as.vector(Y), stress.func, delta=stats::as.dist(Dmat))$par
  Y <- matrix(y, ncol=2)

  plot (Y[,1], Y[,2], asp=1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")

  cluster.vec <- x$clustering
  k <- nlevels(factor(cluster.vec))
  if (missing(col)) col <- grDevices::rainbow(k)
  if (length(col)<k) col <- rep(col,k)

  for (i in 1:k)
    graphics::points(Y[cluster.vec==i,1], Y[cluster.vec==i,2], col=col[i], ...)

  Y
}

# ---------------------------------------------------------------------------------
# boot.silhouette
# ---------------------------------------------------------------------------------

#' Constructs a silhouette plot based on proportion of times items cluster together
#'
#' @param clusboot.out an object of class clusboot
#' @param ... more arguments to be passed to `plot()`
#'
#' @return vector of silhouette widths for each of the clusters
#' @export
#'
#' @examples
#' out <- clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' boot.silhouette(out)

boot.silhouette <- function(clusboot.out, ...)
{
  cluster.vec <- clusboot.out$clustering
  k <- nlevels(factor(cluster.vec))
  Pmat <- clusboot.out$proportions

  sil <- rep(NA,k)
  for (i in 1:k)
  {
    current.clus <- (1:length(cluster.vec))[cluster.vec==i]
    current.p <- Pmat[current.clus, current.clus]
    own.p <- mean(current.p[lower.tri(current.p)])
    other.p <- 0
    for (j in (1:k)[-i])
    {  other.clus <- (1:length(cluster.vec))[cluster.vec==j]
    other.mat <- Pmat[current.clus,other.clus]
    other.mean <- mean(other.mat)
    if (other.mean>other.p) other.p <- other.mean
    }
    sil[i] <- own.p-other.p
  }
  names(sil) <- levels(factor(cluster.vec))
  graphics::barplot (sil, names.arg=1:k, horiz=T, xlim=c(0,1), ...)
  sil
}

# ---------------------------------------------------------------------------------
# complete.linkage
# ---------------------------------------------------------------------------------

#' Wrapper function for performing complete linkage clustering
#'
#' @param X samples x variables data matrix
#' @param k number of clusters
#'
#' @return vector of cluster memberships
#' @export
#'
#' @examples
#' complete.linkage(scale(case.study.psychiatrist), k=6)

complete.linkage <- function (X, k)
{ stats::cutree(stats::hclust(stats::dist(X)), k) }
