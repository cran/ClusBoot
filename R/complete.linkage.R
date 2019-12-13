complete.linkage <-
function (X, k) cutree(hclust(dist(X)), k)
