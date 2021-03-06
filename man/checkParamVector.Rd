% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulatePhyloPoissonLogNormal.R
\name{checkParamVector}
\alias{checkParamVector}
\title{Check Vector Parameter}
\usage{
checkParamVector(x, name, tree)
}
\arguments{
\item{x}{vector of parameters being tested.}

\item{name}{name of the parameter.}

\item{tree}{A phylogenetic tree with n tips.}
}
\description{
Check that the parameters are compatible with the tree. Throws an error if not.
}
\keyword{internal}
