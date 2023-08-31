library(fpc)
set.seed(20000)
face <- rFace(50,dMoNo=2,dNoEy=0,p=2)

test_that("kmeanCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=kmeansCBI, k=5),
                  "clusboot")
})

test_that("additional arguments", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=kmeansCBI, k=5, algorithm="MacQueen"),
                  "clusboot")
})

test_that("distance input", {
  expect_s3_class(clusboot(dist(face), B=50, clustering.func=fpc.clusterboot, clustermethod=kmeansCBI, k=5),
                  "clusboot")
})

test_that("hclustCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=hclustCBI, k=2, method="complete"),
                  "clusboot")
})

test_that("hclusttreeCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=hclusttreeCBI, k=2, method="single"),
                  "clusboot")
})

test_that("disthclustCBI", {
  expect_s3_class(clusboot(dist(face), B=50, clustering.func=fpc.clusterboot, clustermethod=disthclustCBI, k=2, method="ward.D"),
                  "clusboot")
})

test_that("noisemclustCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=noisemclustCBI,
                           k=2, multipleboot=FALSE),
                  "clusboot")
})

test_that("distnoisemclustCBI", {
  expect_s3_class(clusboot(dist(face), B=50, clustering.func=fpc.clusterboot, clustermethod=noisemclustCBI,
                           k=5, multipleboot=FALSE),
                  "clusboot")
})

test_that("claraCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=claraCBI, k=5),
                  "clusboot")
})

test_that("pamkCBI", {
  expect_s3_class(clusboot(dist(face), B=50, clustering.func=fpc.clusterboot, clustermethod=pamkCBI),
                  "clusboot")
})

test_that("dbscanCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=dbscanCBI, eps=1e-5, MinPts=5),
                  "clusboot")
})

test_that("mahalCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=mahalCBI),
                  "clusboot")
})

test_that("mergenormCBI", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, clustermethod=mergenormCBI, G=10, modelNames="EEE",nnk=2),
                  "clusboot")
})

test_that("bootmethod=subset", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, bootmethod="subset", clustermethod=claraCBI, k=5),
                  "clusboot")
})

test_that("bootmethod=subset with dist", {
  expect_s3_class(clusboot(dist(face), B=50, clustering.func=fpc.clusterboot, bootmethod="subset", clustermethod=claraCBI, k=5),
                  "clusboot")
})

test_that("bootmethod=noise", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, bootmethod="noise", clustermethod=claraCBI, k=5),
                  "clusboot")
})

test_that("bootmethod=jitter", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, bootmethod="jitter", clustermethod=claraCBI, k=5),
                  "clusboot")
})

test_that("bootmethod=bojit", {
  expect_s3_class(clusboot(face, B=50, clustering.func=fpc.clusterboot, bootmethod="bojit", clustermethod=claraCBI, k=5),
                  "clusboot")
})
