\name{complete.linkage}
\alias{complete.linkage}
\title{
Wrapper function to return only the clustering allocation.
}
\description{
Wrapper function to return only the clustering allocation from the hclust function with default method = "complete".
}
\usage{
complete.linkage(X, k)
}
\arguments{
  \item{X}{
numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
}
  \item{k}{
either the number of clusters, say k, or a set of initial (distinct) cluster centres. If a number, a random set of (distinct) rows in x is chosen as the initial centres.
}
}
\value{
a single vector of cluster allocations, the output of cutree(hclust(dist(X)), k)
}
\references{
hclust()
}
\author{
Sugnet Lubbe slubbe@sun.ac.za
}
\examples{
data(case_study_psychiatrist)
complete.linkage(scale(case_study_psychiatrist), k=6)
}

