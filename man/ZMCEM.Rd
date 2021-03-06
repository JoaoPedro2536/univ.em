% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_mcem_logll.R
\name{zmcem}
\alias{zmcem}
\title{E-Step for Mu estimation, Simulating Z's}
\usage{
zmcem(theta, data)
}
\arguments{
\item{theta}{The current state of the parameters, vector containing
mean and sd.}

\item{data}{Contingency table, matrix format of three columns, first column
lower limits, second column upper limits, and third column observed
frequencies.}
}
\value{
Simulate data about each specific interval and assign it to the 
 simulate matrix.
}
\description{
E-Step for Mu estimation, Simulating Z's
}
\examples{

 library(truncnorm)

 simdataaaa <- em.univ::univ_simul(ncol_matrix=1,
                         n=50,
                         nclass = 10,
                         mean = 68,
                         sd = 1.80,
                         fr_breaks=c(62,64,66,68,70,72,74,76,78))
 
 zm_matriz <- em.univ::zmcem(theta=c(67,2),data = simdataaaa$simul_data[,,1])
 zm_matriz
}
