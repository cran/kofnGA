\name{kofnGA-package}
\alias{kofnGA-package}
\docType{package}
\title{A genetic algorithm for selection of fixed-size subsets}
\description{
Genetic algorithm to do subset selection:  search for a subset of a fixed size, k, from the 
integers 1:n, such that user-supplied function is minimized at that subset.
}
\details{
\tabular{ll}{
Package: \tab kofnGA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-23\cr
License: \tab GPL-2\cr
}

This package provides the function \code{kofnGA}, which implements a genetic algorithm (GA)
to perform subset selection; that is, choosing the best \emph{k} elements from among \emph{n}
possibilities. We label the set of possibilities from which we are choosing by the integers 
\code{1:n}, and a solution is represented by an index vector, i.e., a vector of integers in 
the range [1, n] (with no duplicates) indicating which members of the set to choose. The 
objective function (defining which solution is \dQuote{best}) is arbitrary and user-supplied; 
the only restriction on this function is that its first argument must be an index vector 
encoding the solution.

The search
results output by \code{kofnGA} are a list object assigned to the S3 class \code{GAsearch}.
The package includes \code{summary}, \code{print}, and \code{plot} methods for this class 
to make it easier to inspect the results.
}
\author{Mark A. Wolters <kramsretlow@gmail.com>}
\references{
Mark A. Wolters (2014, in press), \dQuote{A Genetic Algorithm for Selection of Fixed-Size Subsets}, with Application 
to Design Problems,'' \emph{Journal of Statistical Software}, VOLUME, PAGES
}
\keyword{package}
\examples{
# See the examples in the kofnGA help.
}
