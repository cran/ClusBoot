\name{clusboot}
\alias{clusboot}
\title{
Performs bootstrap on a cluster analysis output}
\description{
B bootstrap samples are drawn with replacement from the data and cluster analysis is performed on the bootstrap samples.
}
\usage{
clusboot(datmat, B = 1000, clustering.func = complete.linkage, ...)
}
\arguments{
  \item{datmat}{
a samples by variables data matrix or a distance object, whatever input is required for clustering.func
}
  \item{B}{
number of bootstrap replicates
}
  \item{clustering.func}{
any function that performs cluster analysis which returns a single vector of cluster allocations as output
}
  \item{\dots}{
any other arguments to be sent to clustering.func
}
}
\value{
  \item{proportions}{an nxn matrix with (i,j)-th element the proportion of times objects i and j clustered in the same cluster)}
  \item{clustering}{vector of cluster allocations of the input data}
}
\references{
Lubbe, S. Visualisations associated with bootstrapping cluster analysis. Data Science, Statistics and Visualisation conference, July 2017, Lisbon, Portugal.
}
\author{
Sugnet Lubbe slubbe@sun.ac.za
}
\seealso{
\code{\link{plot.clusboot}} to visually represent the bootstrap replications and \code{\link{boot.silhouette}} for a silhouette summary of the bootstrap replicates as well as \code{\link{kmeans}}
}
\examples{
data(case_study_psychiatrist)
clusboot (scale(case_study_psychiatrist), B=100, k=6, clustering.func=complete.linkage)
}
