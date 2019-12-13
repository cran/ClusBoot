\name{boot.silhouette}
\alias{boot.silhouette}
\title{
Construct silhouette plot from bootstrap replicates
}
\description{
A silhouette plot is constructed, indicating the proportion of times the cluster members cluster together in the same cluster.
}
\usage{
boot.silhouette(clusboot.out, ...)
}
\arguments{
  \item{clusboot.out}{
an object of class clusboot, usually, a result of a call to clusboot
}
  \item{\dots}{
more plotting parameters, e.g. col
}
}
\details{
The clustering tightness is computed for each cluster, compared to the 'nearest' alternative cluster. The cluster tightness is computed 
as the mean of the proportion of times each pair of objects are clustered together in the same cluster. The tightness to other clusters
is computed as the mean of the proportion an item of this cluster and one from the other cluster as clustered together. The 'nearest' 
alternative cluster is defined as the alternative cluster with the numerically largest tightness to the current cluster.
}
\value{
a vector with number of components equal to the number of classes, providing the silhouette value for each class.
}
\references{
Lubbe, S. Visualisations associated with bootstrapping cluster analysis. Data Science, Statistics and Visualisation conference, July 2017, Lisbon, Portugal.
}
\author{
Sugnet Lubbe slubbe@sun.ac.za
}
\seealso{
\code{\link{clusboot}} for performing bootstrap on a cluster analysis output and \code{\link{plot.clusboot}} to visually represent the bootstrap replications
}
\examples{
data(case_study_psychiatrist)
boot.out <- clusboot (scale(case_study_psychiatrist), B=100, k=6, clustering.func=complete.linkage)
plot(boot.out)
boot.silhouette (boot.out)
}

