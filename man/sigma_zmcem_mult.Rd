% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_bivariado_multivariado.R
\name{sigma_zmcem_mult}
\alias{sigma_zmcem_mult}
\title{Estimate Covariance Matrix for Sample}
\usage{
sigma_zmcem_mult(data, generated_sample, mu_vector)
}
\arguments{
\item{data}{Data in the form of multivariated grouped data.}

\item{generated_sample}{Array, which is the simulated (generated) data.}

\item{mu_vector}{The mean vector from the Mu.}
}
\value{
Return the covariance matrix.
}
\description{
This function calculate the covariance matrix.
}
\examples{


library(MASS)
library(em.univ)
library(tmvtnorm)

simulateddata = em.univ::mult_simul(mm = c(68,68),
                          ss = base::matrix(c(3,2,2,6),2,2) ,
                          n_data_sets = 1,
                          breaks_x = c(-Inf,64,65,66,67,68,69,70,71,72,Inf),
                          breaks_y = c(-Inf,64.2,65.2,66.2,
                                       67.2,68.2,69.2,70.2,
                                       71.2,72.2,Inf),
                          lower_x = base::rep(c(-Inf,64,65,66,67,
                                          68,69,70,71,72),10),
                          lower_y = c(base::rep(-Inf,10),
                                      base::rep(64.2,10),
                                      base::rep(65.2,10),
                                      base::rep(66.2,10),
                                      base::rep(67.2,10),
                                      base::rep(68.2,10),
                                      base::rep(69.2,10),
                                      base::rep(70.2,10),
                                      base::rep(71.2,10),
                                      base::rep(72.2,10)),
                          upper_x = base::rep(c(64,65,66,67,68,69,70,71,72,Inf),10),
                          upper_y = c(base::rep(64.2,10),
                                      base::rep(65.2,10),
                                      base::rep(66.2,10),
                                      base::rep(67.2,10),
                                      base::rep(68.2,10),
                                      base::rep(69.2,10),
                                      base::rep(70.2,10),
                                      base::rep(71.2,10),
                                      base::rep(72,2,10),
                                      base::rep(Inf,10))
)

mu <- c(67,67)
sigma <- base::matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)

generated_sample<- em.univ::zmcem_mult(data=simulateddata[,,1],
                              mu=mu,
                              sigma=sigma)

mu_vector<- em.univ::mu_zmcem_mult(data=simulateddata[,,1],
                          generated_sample=generated_sample)


out<- em.univ::sigma_zmcem_mult(data=simulateddata[,,1],
                       generated_sample = generated_sample,
                       mu_vector = mu_vector)

out

}
