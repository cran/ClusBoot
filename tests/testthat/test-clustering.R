X <- iris[,1:4]
Dmat <- dist(X)

test_that("complete linkage - data matrix", {
  expect_equal(cutree(hclust(dist(X)),k=3),
               complete.linkage(X,k=3))
})

test_that("clustering results are the same - data matrix", {
  expect_equal(cutree(hclust(dist(X)),k=3),
               clusboot(X,B=2,clustering.func=complete.linkage,k=3)$clustering)
})

complete.dist.wrapper <- function(Dobj, k)
  { cutree(hclust(Dobj), k) }

test_that("clustering results are the same - dist matrix", {
  expect_equal(cutree(hclust(dist(iris[,1:4])),k=3),
               clusboot(Dmat,B=2,clustering.func=complete.dist.wrapper,k=3)$clustering)
})

