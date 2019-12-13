boot.silhouette <-
function(clusboot.out, ...)
{
  cluster.vec <- clusboot.out[[2]]
  k <- nlevels(factor(cluster.vec))
  Pmat <- clusboot.out[[1]]

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
  barplot (sil, names.arg=1:k, horiz=T, xlim=c(0,1), ...)
  sil
}
