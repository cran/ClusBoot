# Sugnet Lubbe
# MuViSU (Centre for Multi-Dimensional Data Visualisation)
# Department of Statistics and Actuarial Science
# Stellenbosch University, June 2022
# ============================================================

# ---------------------------------------------------------------------------------
# fpc.clusterboot
# function adapted from fpc package to perform different bootstrapping options
# ---------------------------------------------------------------------------------
#' Resampling according to the methods discussed in Hennig (2007)
#'
#' @param data a data matrix or distance object which will be the input to the clustering function
#' @param B number of bootstrap replicates
#' @param distances see `?fpc::clusterboot`
#' @param bootmethod see `?fpc::clusterboot`
#' @param bscompare see `?fpc::clusterboot`
#' @param multipleboot see `?fpc::clusterboot`
#' @param jittertuning see `?fpc::clusterboot`
#' @param noisetuning see `?fpc::clusterboot`
#' @param subtuning see `?fpc::clusterboot`
#' @param clustermethod see `?fpc::clusterboot`
#' @param noisemethod see `?fpc::clusterboot`
#' @param count see `?fpc::clusterboot`
#' @param seed see `?fpc::clusterboot`
#' @param datatomatrix see `?fpc::clusterboot`
#' @param ... additional arguments to be sent to the function specified in clustermethod
#'
#' @export
#' @return a list with two components; boot.out contains the computations for clusboot and out contains the clustering
#'           solution of the original data set
#'
fpc.clusterboot <- function(data, B, distances = (inherits(data, "dist")),
                            bootmethod = "boot", bscompare = TRUE, multipleboot = FALSE,
                            jittertuning = 0.05, noisetuning = c(0.05, 4), subtuning = floor(nrow(data)/2),
                            clustermethod, noisemethod = FALSE, count = TRUE, seed = NULL, datatomatrix = TRUE,
                            ...)
{
      sumlogic <- function(x, y, relation = "eq") switch(relation,
                                                       eq = sum(x == y, na.rm = TRUE), s = sum(x < y, na.rm = TRUE),
                                                       l = sum(x > y, na.rm = TRUE), se = sum(x <= y, na.rm = TRUE),
                                                       le = sum(x >= y, na.rm = TRUE))
   if (!is.null(seed))
      set.seed(seed)
    invisible(distances)
    if (datatomatrix)
      data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    if (datatomatrix & !distances) {
      cod <- stats::cov(data)
      md <- colMeans(data)
    }
    lb <- length(bootmethod)

    if (distances)
      c1 <- clustermethod(stats::as.dist(data), ...)
    else c1 <- clustermethod(data, ...)
    if (noisemethod) {
      if (c1$nccl == 0)
        stop("No clusters, only noise estimated!")
    }
    else c1$nccl <- c1$nc

    if (("jitter" %in% bootmethod) | ("bojit" %in% bootmethod)) {
      if (!datatomatrix | distances)
        stop("datatomatrix=FALSE and distances require boot or subset as bootmethod.")
      jsd <- numeric(0)
      ecd <- eigen(cod, symmetric = TRUE)
      ecd$values[ecd$values < 0] <- 0
      ecd$values[is.na(ecd$values)] <- 0
      rotdata <- data %*% solve(t(ecd$vectors))
      for (i in 1:p) {
        sx <- sort(rotdata[, i])
        dx <- sx[2:n] - sx[1:(n - 1)]
        dx <- dx[dx > 0]
        jsd[i] <- stats::quantile(dx, jittertuning)
      }
    }
    if ("noise" %in% bootmethod) {
      if (!datatomatrix | distances)
        stop("datatomatrix=FALSE and distances require boot or subset as bootmethod.")
      ecd <- eigen(cod, symmetric = TRUE)
      ecd$values[ecd$values < 0] <- 0
    }

    boot.out <- matrix (NA, nrow=2*nrow(data)*nrow(data), ncol=B)
# *************************************************
      for (b in 1:B) {
        if (count) cat(bootmethod[1], b, "\n")
        if (bootmethod[1] == "boot")
            {  bsamp <- sample(n, n, replace = TRUE)
               if (!multipleboot) bsamp <- unique(bsamp)
               if (distances) mdata <- data[bsamp, bsamp]
               else mdata <- data[bsamp, ]
            }
        if (bootmethod[1] == "subset")
            {  bsamp <- sample(n, subtuning, replace = FALSE)
               if (distances) mdata <- data[bsamp, bsamp]
               else mdata <- data[bsamp, ]
            }
        if (bootmethod[1] == "jitter")
            {  jnoise <- matrix(0, ncol = p, nrow = n)
               for (j in 1:p) jnoise[, j] <- stats::rnorm(n, sd = jsd[j])
               jnoise <- jnoise %*% t(ecd$vectors)
               mdata <- data + jnoise
               bsamp <- 1:n
            }
        if (bootmethod[1] == "bojit")
            {  bsamp <- sample(n, n, replace = TRUE)
               jnoise <- matrix(0, ncol = p, nrow = n)
               for (j in 1:p) jnoise[, j] <- stats::rnorm(n, sd = jsd[j])
               jnoise <- jnoise %*% t(ecd$vectors)
               mdata <- data[bsamp, ] + jnoise
            }
        if (bootmethod[1] == "noise")
            {  noiseind <- as.logical(stats::rbinom(n, 1, noisetuning[1]))
               nn <- sum(noiseind)
               jnoise <- matrix(0, ncol = p, nrow = nn)
               for (j in 1:p) jnoise[, j] <- stats::runif(nn, min = -noisetuning[2] *
                                                sqrt(ecd$values[j]), max = noisetuning[2] *
                                                sqrt(ecd$values[j]))
               jnoise <- t(t(jnoise %*% t(ecd$vectors)) + md)
               mdata <- data
               mdata[noiseind, ] <- jnoise
               bsamp <- (1:n)[!noiseind]
             }

        if ("diss" %in% names(formals(clustermethod)) & distances)
          bc1 <- clustermethod(mdata, diss = TRUE, ...)
        else bc1 <- clustermethod(mdata, ...)
        if (noisemethod) {
          effnc1 <- c1$nccl
          effnb1 <- bc1$nccl
        }
        else {
          effnc1 <- c1$nc
          effnb1 <- bc1$nc
        }

        if ("diss" %in% names(formals(clustermethod)) & distances)
          bc1 <- clustermethod(mdata, diss = TRUE, ...)
        else bc1 <- clustermethod(mdata, ...)

        out <- rep(0,n)
        for (i in 1:effnb1)
          out[bc1$clusterlist[[i]]] <- i
        clusD <- totD <- matrix (0, nrow=n, ncol=n)

        boot.sample <- table(bsamp)
        boot.names <- as.numeric(names(boot.sample))
        Dmat <- matrix (boot.sample,ncol=1) %*% matrix (boot.sample, nrow=1)
        totD[boot.names,boot.names] <- Dmat

        kk <- nlevels(factor(out))
          for (i in 1:kk)
          {
            clus.count <- table(bsamp[out==i])
            clus.names <- as.numeric(names(clus.count))
            Dmat <- matrix (clus.count,ncol=1) %*% matrix (clus.count,nrow=1)
            clusD[clus.names,clus.names] <- Dmat
          }

        boot.out[,b] <- c(clusD, totD)
      }
    # *************************************************

    out <- rep(0,n)
    for (i in 1:length(c1$clusterlist))
      out[c1$clusterlist[[i]]] <- i
    list(boot.out=boot.out, out=out)
}

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
#' @return an object of class `clusboot` which is a list with the following components:
#' \item{proportions}{matrix of size nxn with cell ij containing the proportion of bootstrap replicates
#'                    in which object i and object j clustered together.}
#' \item{clustering}{a vector of length n containing the cluster membership of the n input objects.}
#' \item{sil}{a vector of length the number of clusters containing the bootstrap-silhouette values for the clusters.}
#' \item{indv.sil}{a vector of length n containing the bootstrap-silhouette values for the individual objects.}
#' \item{sil.order}{a vector of length n containing the ordering of the n objects used by the functions
#'                  `boot.silhouette` and `boot.proportions` to order objects in the same cluster adjacent and
#'                  clusters in decreasing order of cluster tightness.}
#' \item{ave.sil.width}{the overall stability of the clustering solution, obtained by averaging over the individual object bootstrap-silhouette values.}
#'
#' @export
#'
#' @details
#' Any R function performing cluster analysis can be specified in `clustering.func` although a wrapper function is
#' typically needed to isolate only the vector output of cluster memberships. See `?complete.linkage` as an example.
#' Should users perfer to use alternative resamling schemes, other than the bootstrap, Hennig (2007) discuss a variety
#' of options which could be accessed by specifying `clustering.func = fpc.clusterboot`. In addition, the sampling
#' method is specified in the argument `bootmethod` and additional arguments for the function `clusterboot` in the
#' package `fpc` must be given. Note that only the resampling facilities of `clusterboot` is utilised while the
#' computation of proportions and silhouette widths remain unchanged. The output object of class `clusboot`
#' will remain unchanged as only the resampling section of `clusterboot` is used.
#'
#' @examples
#' clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' library(fpc)
#' clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=fpc.clusterboot,
#'           clustermethod=hclustCBI, method="complete", bootmethod="subset", subtuning=10)
#'
#' @references Hennig, C., 2007. Cluster-wise assessment of cluster stability. Computational Statistics & Data Analysis,
#' 52(1), pp.258-271.

clusboot <- function(datmat, B=1000, clustering.func=complete.linkage, ...)
{
  obj <- datmat
  if (data.class(datmat) != "matrix") datmat <- as.matrix(datmat)
  n <- nrow(datmat)
  if (is.null(rownames(datmat))) sample.names <- 1:n
  else sample.names <- rownames(datmat)
  boot.samples <- matrix(sample (1:n, size=n*B, replace=T), ncol=B)

  if (!identical(clustering.func, fpc.clusterboot))
  {
  boot.out <- apply(boot.samples,2,function(x)
                                   { if (data.class(obj) == "dist")
                                          out <- clustering.func(stats::as.dist(datmat[x,x]), ...)
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
  }
  else
  {
    fpc.boot.out <- fpc.clusterboot(data=obj, B = B, ...)
    boot.out <- fpc.boot.out$boot.out
  }

  boot.out <- apply (boot.out, 1, sum)
  clusD <- matrix (boot.out[1:(n*n)], nrow=n, ncol=n)
  totD <- matrix (boot.out[-(1:(n*n))], nrow=n, ncol=n)

  clusD <- clusD/totD
  dimnames(clusD) <- list(sample.names,sample.names)
  if (!identical(clustering.func, fpc.clusterboot))
    out <- clustering.func(obj, ...)
  else out <- fpc.boot.out$out
  output <- list(proportions=clusD, clustering=out)
  class(output) <- "clusboot"
  output <- calc.silhouette(output)
  output
}

# ---------------------------------------------------------------------------------
# plot.clusboot
# ---------------------------------------------------------------------------------

#' MDS plot of similarities given by the proportion of bootstrap replicates where objects cluster together
#'
#' @param x an object of class clusboot
#' @param col single colour or a vector specifying a colour for each object
#' @param show.silhouette logical indicating whether plotting character size should represent the inidivual silhouette values
#' @param ... more arguments to be passed to `plot()`
#'
#' @return matrix of similarities (proportions)
#' @export
#'
#' @examples
#' out <- clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' plot(out)

plot.clusboot <- function(x, col, show.silhouette=TRUE, ...)
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
    if (show.silhouette) graphics::points(Y[cluster.vec==i,1], Y[cluster.vec==i,2], col=col[i], cex=(1-x$indv.sil[cluster.vec==i])*3+0.1, ...)
    else graphics::points(Y[cluster.vec==i,1], Y[cluster.vec==i,2], col=col[i], ...)
  Y
}

# ---------------------------------------------------------------------------------
# boot.proportions
# ---------------------------------------------------------------------------------

#' Heatmap of the proportion of bootstrap replicates where objects cluster together
#'
#' @param x an object of class clusboot
#' @param col vector of colours for shading to indicate proportion values
#' @param show.vals logical value indicating whether proportion values should be added to individual cells
#' @param text.col colour of text for show.vals if `TRUE`
#' @param cluster.col colour of lines demarcating cluster membership
#' @param ... more arguments to be passed to `plot()`
#'
#' @export
#'
#' @examples
#' out <- clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' boot.proportions(out)

boot.proportions <- function(x, col=grDevices::heat.colors(101,rev=TRUE), show.vals=F, text.col="black", cluster.col="firebrick", ...)
{
  Pmat <- x$proportions[x$sil.order, x$sil.order]
  clus.vec <- x$clustering[x$sil.order]
  n <- nrow(Pmat)
  m <- length(col)-1
  old.par <- graphics::par(mar=c(2.1, 2.1, 0.1, 0.1), no.readonly=TRUE)
  on.exit(graphics::par(old.par))
  plot (0:n, 0:n, type="n", xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", ...)
  graphics::axis (side=1, at=(1:n)-0.5, label=colnames(Pmat), las=2, cex.axis=0.6)
  graphics::axis (side=2, at=(n:1)-0.5, label=rownames(Pmat), las=2, cex.axis=0.6)
  for (x in 1:n)
    for (y in 1:n)
      { graphics::rect(x-1, n-y+1, x, n-y, col=col[floor(Pmat[x,y]*m+1)], border=NA)
        if (show.vals) graphics::text (x-0.5, n-y+0.5, format(round(Pmat[x,y]*100)/100,nsmall=2), cex=0.4, adj=c(0.5,0.5))
    }
  for (i in 2:n)
    if (clus.vec[i-1] != clus.vec[i])
      {
         graphics::lines (rep(i-1,2), c(0,n), col=cluster.col)
         graphics::lines (c(0,n), rep(n-i+1,2), col=cluster.col)
      }
  invisible(NULL)
}

# ---------------------------------------------------------------------------------
# calc.silhouette
# ---------------------------------------------------------------------------------

#' Computes the silhouette vales based on proportion of times items cluster together
#'
#' @param clusboot.out an object of class clusboot
#'
#' @return an object of class clusboot

calc.silhouette <- function(clusboot.out)
{
  cluster.vec <- clusboot.out$clustering
  k <- nlevels(factor(cluster.vec))
  cluster.names <- levels(factor(cluster.vec))
  n <- length(cluster.vec)
  Pmat <- clusboot.out$proportions

  sil <- rep(NA,k)
  indv.sil <- rep(NA, n)
  for (i in 1:k)
  {
    current.clus <- (1:length(cluster.vec))[cluster.vec==cluster.names[i]]
    current.p <- Pmat[current.clus, current.clus,drop=F]
    temp.p <- current.p
    diag(temp.p) <- NA
    indv.sil[current.clus] <- apply(temp.p, 2, mean, na.rm=T)
    own.p <- mean(current.p[lower.tri(current.p)])
    other.p <- 0
    for (j in (1:k)[-i])
      {  other.clus <- (1:length(cluster.vec))[cluster.vec==j]
         other.mat <- Pmat[current.clus,other.clus]
         other.mean <- mean(other.mat)
         if (!is.na(other.mean)) if (other.mean>other.p) other.p <- other.mean
      }
    sil[i] <- own.p-other.p
    indv.sil[current.clus] <- indv.sil[current.clus] - other.p
  }
  names(sil) <- cluster.names
  names(indv.sil) <- paste(cluster.vec,rownames(Pmat),sep="-")

  sil.order <- cluster.names[order(sil)]
  indv.sil.order <- order(match(cluster.vec,sil.order), indv.sil, decreasing=TRUE)
  clusboot.out$sil <- sil
  clusboot.out$indv.sil <- indv.sil
  clusboot.out$sil.order <- indv.sil.order
  clusboot.out$avg.sil.width <- mean(indv.sil,na.rm=TRUE)
  clusboot.out
}

# ---------------------------------------------------------------------------------
# boot.silhouette
# ---------------------------------------------------------------------------------

#' Produces silhouette plots
#'
#' @param clusboot.out an object of class clusboot
#' @param ... more arguments to be passed to `barplot()`
#'
#' @export
#'
#' @return list of silhouette widths
#'
#' @examples
#' out <- clusboot (scale(case.study.psychiatrist), B=100, k=6, clustering.func=complete.linkage)
#' boot.silhouette(out)

boot.silhouette <- function(clusboot.out, ...)
{
  old.par <- graphics::par(mfrow=c(1,2), no.readonly = TRUE)
  on.exit (graphics::par(old.par))
  sil.val <- sort(clusboot.out$sil, na.last=TRUE)
  indv.val <- clusboot.out$indv.sil
  n <- length(indv.val)
  clus.num <- clusboot.out$clustering[rev(clusboot.out$sil.order)]

  mp <- graphics::barplot(indv.val[rev(clusboot.out$sil.order)], horiz=T, las=2, cex.names=0.5, ...)
  for (i in 2:length(mp))
    if (clus.num[i-1] != clus.num[i]) graphics::lines (graphics::par("usr")[1:2], rep((mp[i-1]+mp[i])/2,2), col="firebrick")
  graphics::barplot(sil.val, horiz=T, xlim=c(0,1), ...)
  list (individual.silhouette.values = clusboot.out$indv.sil, cluster.silhouette.values = clusboot.out$sil)
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
