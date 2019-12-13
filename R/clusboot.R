clusboot <-
function (datmat, B=1000, clustering.func=complete.linkage, ...)
{
  n <- nrow(datmat)
  if (is.null(rownames(datmat))) rownames(datmat) <- 1:n
  boot.samples <- matrix(sample (1:n, size=n*B, replace=T), ncol=B)

  boot.out <- apply(boot.samples,2,function(x) 
                                     { out <- clustering.func(datmat[x,], ...)

                                       clusD <- totD <- matrix (0, nrow=n, ncol=n, dimnames=list(rownames(datmat),rownames(datmat)))

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
  dimnames(clusD) <- list(rownames(datmat),rownames(datmat))
  out <- clustering.func(datmat, ...)
  output <- list(proportions=clusD[order(out),order(out)], clustering=out[order(out)])
  class(output) <- "clusboot"
  output
}
