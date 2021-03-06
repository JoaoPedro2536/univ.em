% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_mcem_logll.R
\name{mcem}
\alias{mcem}
\title{This is the maximization step of the mcem algorithm (M-step).}
\usage{
mcem(data, theta_init, maxit = 1000, tol1 = 0.01, tol2 = 0.001)
}
\arguments{
\item{data}{Contingency table, matrix format of three columns, first column
lower limits, second column upper limits, and third column observed
frequencies.}

\item{theta_init}{The initial values of the parameters, for mu and sigma. 
Vector with two values.}

\item{maxit}{The maximum number of iteration of the EM algorithm.}

\item{tol1}{A number, the stopping criteria for updating mu.}

\item{tol2}{A number, the stopping criteria for updating sigma.}
}
\value{
This is the maximization step of the mcem algorithm (M-step) that has
 defined it using the function E-Step for mean estimate and variance
 estimate. Returns the estimates for mean (mu) and sigma (variance).
}
\description{
This is the maximization step of the mcem algorithm (M-step).
}
\examples{

 library(em.univ)


 simdataaaa <- em.univ::univ_simul(ncol_matrix=1,
                         n=50,
                         nclass = 10,
                         mean = 68,
                         sd = 1.80,
                         fr_breaks=c(62,64,66,68,70,72,74,76,78))
 
 outputmcem2 <- base::list()
 
 outputmcem2[[1]] <- em.univ::mcem(data=simdataaaa$simul_data[,,1],
                           theta_init=c(67,2),
                           maxit = 1000,tol1=1e-2,tol2=1e-3)
 
 outputmcem2
}
