% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NormalBiVarMod}
\alias{NormalBiVarMod}
\title{Modified Simulation of two groups of normal bivariate distributions}
\format{A data frame 350 rows and 3 variables: \describe{ \item{Y}{factor,
  response variable used to codify the groups, also known as labels}
  \item{X1}{numeric, first covariable used for clasification}
  \item{X1}{numeric, second covariable used for clasification}}}
\usage{
NormalBiVarMod
}
\description{
A dataset containing two groups of normal bivariate distributions with
diferent means and covariance matriz used to test the bpwpm2 model. Group 1
has 200 observations with mean (2,2) and covariance matrix (0.25, 0.35; 0.35,
1). Group 2, 150 observations, with mean (4,4) and covariance matrix (1,
-0.24; -0.24, 0.64). Afterwards, the 13th right-most observation of Group 2
were manually re-coded to belong to Group 1 in order to have disjoined groups
to make the classification problem more intresting.
}
\keyword{datasets}
