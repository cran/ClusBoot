\name{plot.clusboot}
\alias{plot.clusboot}
\title{
Constructs a Multi-dimensional scaling plot to represent the bootstrap replications.
}
\description{
Constructs a Multi-dimensional scaling plot to represent the bootstrap replications.
}
\usage{
\method{plot}{clusboot}(x, col = NULL, ...)
}
\arguments{
  \item{x}{
an object of class clusboot, usually, a result of a call to clusboot
}
  \item{col}{
an optional vector equal in length to the number of clusters with colours to represent different clusters
}
  \item{\dots}{
more plotting parameters
}
}
\details{
Constructs a Multi-dimensional scaling plot to represent the bootstrap replications.
}
\value{
a two-column matrix with rows containing the MDS coordinates of the samples 
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
}
