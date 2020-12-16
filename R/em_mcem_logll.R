
### E-Step For Mean Estimate ---------------------------------------------------

#' E-Step For Mean Estimate
#'
#' @param theta Initial theta value, vector containing mean and sd.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param Freq Frequency over the intervals, values in vector.
#' @return Return M which are the estimates of mean in E-step
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  Munew <- Mest(theta=c(67,2),
#'                bl=simdataaaa$simul_data[,1,1],
#'                bu=simdataaaa$simul_data[,2,1],
#'                Freq=simdataaaa$simul_data[,3,1])
#'  
#'  
#'  Munew



Mest<- function(theta,bl,bu,Freq){

  Aj<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))

  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]
  }

  for(i in 1:base::length(bl)){
    Aj[i]<- theta[1] - theta[2]*
      ((stats::dnorm(bstar[i]) - stats::dnorm(astar[i]))/
         ((stats::pnorm(bstar[i]) - stats::pnorm(astar[i]))) )
  }

  M <- base::sum(Aj*Freq)/sum(Freq)
  return(M)

}
################################################################################

### E-Step For Variance Estimate -----------------------------------------------

#' E-Step For Variance Estimate
#'
#' @param theta Initial theta value, vector containing mean and sd.
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param Freq Frequency over the intervals, values in vector.
#' @param muupdate Is the updated estimates of mu.
#' @return Return SS which are the estimates of sigma (variance) in E-step.
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                         n=50,
#'                         nclass = 10,
#'                         mean = 68,
#'                         sd = 1.80,
#'                         fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  SSnew<- SSest(theta=c(67,2),
#'               bl=simdataaaa$simul_data[,1,1],
#'               bu=simdataaaa$simul_data[,2,1],
#'               muupdate=Mest(theta=c(67,2),
#'                             bl=simdataaaa$simul_data[,1,1],
#'                             bu=simdataaaa$simul_data[,2,1],
#'                             Freq=simdataaaa$simul_data[,3,1]),
#'               Freq=simdataaaa$simul_data[,3,1])
#'  
#'  SSnew


SSest<- function(theta,bl,bu,muupdate,Freq){

  Bj<- base::rep(0,base::length(bl))
  bstar<- base::rep(0,base::length(bl))
  astar<- base::rep(0,base::length(bl))

  for(i in 1:base::length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2]
    astar[i]<- (bl[i]-theta[1])/theta[2]

  }

  astar[1] <- (-1000)
  bstar[base::length(bl)]<- (1000)


  for(i in 1:base::length(bl)){
    Bj[i]<- theta[2]^2*
      (1-(bstar[i]*stats::dnorm(bstar[i]) - astar[i]*stats::dnorm(astar[i]))/
         (stats::pnorm(bstar[i]) - stats::pnorm(astar[i])))+
      (muupdate - theta[1])^2 + (2*theta[2]*(muupdate-theta[1])*
         ((stats::dnorm(bstar[i]) - stats::dnorm(astar[i]))/
            (stats::pnorm(bstar[i]) - stats::pnorm(astar[i]))))

  }

  SS<- base::sum(Bj*Freq)/base::sum(Freq)
  return(SS)
}
################################################################################

### M-step maximization step of the EM algorithm -------------------------------

#' Maximization step of the EM algorithm
#'
#' @param bl Lower bound of the intervals, values in vector, starting from -inf.
#' @param bu Upper bound of the intervals, values in vector, ending with +inf.
#' @param Freq Frequency over the intervals, values in vector.
#' @param theta_init The initial values of the parameters, for mu and sigma. 
#'  Vector with two values.
#' @param maxit The maximum number of iteration of the EM algorithm.
#' @param tol1 A number, the stopping criteria for updating mu.
#' @param tol2 A number, the stopping criteria for updating sigma.
#' @return This is the maximization step of the EM algorithm (M-step) that has
#'  defined it using the function E-Step for mean estimate and variance
#'  estimate. Return a list has as arguments "mu_estimate" for the average and
#'  "sigma_estimate" for the variance.
#' @export
#' @examples
#'  output2 =list()
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                        n=50,
#'                        nclass = 10,
#'                        mean = 68,
#'                        sd = 1.80,
#'                        fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'
#'  output2<- EM(bl=simdataaaa$simul_data[,1,1],
#'                bu=simdataaaa$simul_data[,2,1],
#'                Freq=simdataaaa$simul_data[,3,1],
#'                theta_init=c(67,2),
#'                maxit = 1000,
#'                tol1=1e-3,
#'                tol2=1e-4)
#'  output2


EM<- function(bl,bu,Freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]

  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    Munew<- Mest(theta=cur,bl,bu,Freq)
    SSnew<- SSest(theta=cur,bl,bu,
                  muupdate=Mest(theta=cur,bl,bu,Freq) ,Freq)

    Mu_new <- Munew
    S_new <- base::sqrt(SSnew)
    new_step<- c(Mu_new,S_new)

    if(base::abs(cur[1]-new_step[1])<tol1 & base::abs(cur[2]-new_step[2])<tol2){
      flag<-1 ;break
      }
    Mu_cur<- Mu_new
    S_cur<- S_new
  }

  if(!flag) warning
  updateres<- base::list("mu_estimate" = Mu_cur,
                         "sigma_estimate" = (S_cur)^2)

  return(updateres)
}

################################################################################

### func univ_Simul (gerando simdata2) -----------------------------------------

#' Generating data from a gaussian
#'
#' @param ncol_matrix A number, number of repetitions of the random simulated
#'  samples.
#' @param n A number, number of observations in each sample.
#' @param nclass A number, number of classes on the contingency table.
#' @param mean A number, parameter mu value to generate data.
#' @param sd A number, parameter sigma value to generate data.
#' @param fr_breaks A vector, vector containing the values of the referenced
#'  as counts.
#' @return Returns a list with two objects, the first "simul_data" returns
#'  a matrix in which the first column the lower limits of the intervals,
#'  in the second column the values of the upper limits of the intervals and
#'  the third column the counts of the observations. The second object "med"
#'  returns a matrix, in which the first column shows the sampled value and
#'  the second column the observation count.
#' @export
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                        n=50,
#'                        nclass = 10,
#'                        mean = 68,
#'                        sd = 1.80,
#'                        fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  simdataaaa


univ_Simul <- function(ncol_matrix=30,
                       n=50,
                       nclass = 10,
                       mean = 68,
                       sd = 1.80,
                       fr_breaks=c(62,64,66,68,70,72,74,76,78)){

  sim2<- base::matrix(rep(0,n*ncol_matrix),
                ncol=ncol_matrix)
  for(i in 1:base::ncol(sim2)){
    sim2[,i]<- stats::rnorm(n = n,
                            mean = mean,
                            sd = sd)
  }

  ### ### ### ###
  Fr<- base::matrix(rep(0,10*ncol(sim2)),
              ncol=ncol(sim2))

  for(i in 1:base::ncol(sim2)){
    Fr[,i]<- base::table(cut(sim2[,i],
                       breaks=c(-Inf,fr_breaks,Inf)))

  }
  ### ### ### ###

  simdata2<- base::array(rep(0,10*3*ncol_matrix),
                   c(nclass,3,ncol_matrix))

  med2<- base::array(rep(0,10*2*ncol_matrix),
               c(nclass,2,ncol_matrix))

  for(i in 1:base::ncol(sim2)){
    simdata2[,1,i]<- c(-Inf,fr_breaks)
    simdata2[,2,i]<- c(fr_breaks,Inf)
    simdata2[,3,i]<- Fr[,i]

    med2[,1,i]<- (simdata2[,1,i]+simdata2[,2,i])/2
    med2[1,1,i]<- base::min(fr_breaks)-1
    med2[nclass,1,i]<- base::max(fr_breaks)+1
    med2[,2,i]<- Fr[,i]
  }

  final_list <- base::list("simul_data" = simdata2,
                    "med" = med2)

  return(final_list)

}
################################################################################





### E-Step for Mu estimation: Simulating Z's  ----------------------------------
###Calculation of updates for Mu###

#' E-Step for Mu estimation, Simulating Z's
#'
#' @param theta The current state of the parameters, vector containing
#'  mean and sd.
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies.
#' @return Simulate data about each specific interval and assign it to the 
#'  simulate matrix.
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  ZM_matriz = ZMCEM(theta=c(67,2),data = simdataaaa$simul_data[,,1])
#'  ZM_matriz

ZMCEM<- function(theta,data){

  k<- 1000
  sim<- base::matrix(rep(0,k*base::nrow(data)),ncol=base::nrow(data))

  for(i in 1:base::nrow(data)) {
    sim[,i]<- truncnorm::rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2])
  }

  return(sim)
}

################################################################################

### E-Step for Mu --------------------------------------------------------------

#' E-Step estimate the parameter of mu.
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies. The current state of the parameters, vector containing
#'  mean and sd.
#' @param simZ The matrix, which is the samples simulated over the intervals,
#'  of the ZMCEM function.
#' @return Return mu which are the estimates of mean in E-step.
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  ZM_matriz = ZMCEM(theta=c(67,2),data = simdataaaa$simul_data[,,1])
#'  
#'  mu_estimate = MuMCEM(data = simdataaaa$simul_data[,,1], simZ = ZM_matriz)
#'  mu_estimate


MuMCEM<- function(data,simZ){

  n<- base::sum(data[,3])
  Z<- base::colMeans(simZ)
  numerator<- base::rep(0,base::nrow(data))
  for (i in 1:base::nrow(data)) {
    numerator[i]<- data[i,3]*Z[i]
  }

  MuN<- (1/n)*base::sum(numerator)
  return(MuN)
}

################################################################################

### Calculate estimate of sigma ------------------------------------------------

#' E-Step estimate the parameter of sigma.
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies. The current state of the parameters, vector containing
#'  mean and sd.
#' @param simZ The matrix, which is the samples simulated over the intervals,
#'  of the ZMCEM function.
#' @param mupd Is the updated estimates of mu.
#' @return Return sigma which are the estimates of variance in E-step.
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  ZM_matriz = ZMCEM(theta=c(67,2),data = simdataaaa$simul_data[,,1])
#'  
#'  mu_estimate = MuMCEM(data = simdataaaa$simul_data[,,1], simZ = ZM_matriz)
#'  
#'  sigma_estimate = sigmaMCEM(data = simdataaaa$simul_data[,,1],
#'                             simZ = ZM_matriz,
#'                             mupd = mu_estimate)
#'  sigma_estimate



sigmaMCEM<- function(data,simZZ,mupd){

  n<- base::sum(data[,3])
  ZZ<- simZZ
  NewZ<- (ZZ-mupd)^2
  SZNEW<- base::colMeans(NewZ)
  numerator<- base::rep(0,nrow(data))

  for (i in 1:base::nrow(data)){
    numerator[i]<- data[i,3]*SZNEW[i]
  }

  sigmaNN<- (1/n)*base::sum(numerator)
  sig<- base::sqrt(sigmaNN)

  return(sig)
}

################################################################################

### MONTE CARLO EM -------------------------------------------------------------

#' This is the maximization step of the MCEM algorithm (M-step).
#'
#' @param data Contingency table, matrix format of three columns, first column
#'  lower limits, second column upper limits, and third column observed
#'  frequencies.
#' @param theta_init The initial values of the parameters, for mu and sigma. 
#'  Vector with two values.
#' @param maxit The maximum number of iteration of the EM algorithm.
#' @param tol1 A number, the stopping criteria for updating mu.
#' @param tol2 A number, the stopping criteria for updating sigma.
#' @return This is the maximization step of the MCEM algorithm (M-step) that has
#'  defined it using the function E-Step for mean estimate and variance
#'  estimate. Returns the estimates for mean (mu) and sigma (variance).
#' @export
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  outputMCEM2 = list()
#'  
#'  outputMCEM2[[1]]<- MCEM(data=simdataaaa$simul_data[,,1],theta_init=c(67,2),
#'                          maxit = 1000,tol1=1e-2,tol2=1e-3)
#'  
#'  outputMCEM2


MCEM<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){
  
  flag<- 0
  Mu_cur<- theta_init[1]
  S_cur<- theta_init[2]
  iter<- base::rep(0,maxit)
  Svec<- base::rep(0,maxit)
  Mvec<- base::rep(0,maxit)
  mydat<- data
  
  for (i in 1:maxit){
    cur<- c(Mu_cur,S_cur)
    Munew<- MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat))

    Snew<- sigmaMCEM(data=mydat,simZZ=ZMCEM(theta=cur,data=mydat),
                     mupd=MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat)))

    Mu_new<- Munew
    S_new<- Snew
    new_step<- c(Mu_new,S_new)

    if(base::abs(cur[1]-new_step[1])<tol1 & base::abs(cur[2]-new_step[2])<tol2){
      flag<-1 ;break
    }

    Mu_cur<- Mu_new
    S_cur<- S_new
    iter[i]<- i
    Svec[i]<- S_new
    Mvec[i]<- Mu_new

  }

  if(!flag) warning("Didn't Converge \n")
  update <- base::list("mu_estimate" = Mu_cur,
                       "sigma_estimate" = (S_cur)^2)

  return(update)
}
################################################################################





### LOG L FUNCTION -------------------------------------------------------------

#' Likelihood Estimation.
#'
#' @param TL Lower end of the intervals.
#' @param freq Frequency on each interval.
#' @param theta The arguments including mu and sigma.
#' @return Return the result of the log likelihood for the mean (mu) and sigma 
#'  (standard deviation) estimates.
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  
#'  TL<- simdataaaa$simul_data[,1,1]
#'  freq<- simdataaaa$simul_data[,3,1]
#'  res <- stats::optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,
#'  method="L-BFGS-B",lower=c(0.02,0.2),upper=c(180,2))
#'  
#'  estimate<- res$par
#'  mu_estimate<- base::round(res$par[1]/res$par[2],4)
#'  sigma_estimate<- base::round(1/res$par[2],4) 
#'  sigma_squared_estimate = sigma_estimate^2
#'  
#'  mu_estimate
#'  sigma_squared_estimate


Logll <- function(TL,freq,theta){
  m<- length(TL)

  if( (pnorm(theta[2]*TL[2]-theta[1])) < 1e-16  ){
    a <- -1e+6
  } #end if

  else{
    a <- freq[1]*log(pnorm(theta[2]*TL[2]- theta[1]))
  } #end else

  if( (1-pnorm(theta[2]*TL[m]-theta[1])) < 1e-16  ) {
    b <- -1e+6
  }#end if

  else{
    b <- freq[m]*log(1-pnorm(theta[2]*TL[m]-theta[1]))
  }#end else

  c<-0
  for(i in 2:(m-1)){

    if ( (pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1])) < 1e-16 ){
      c <- c -1e+6
    }#end if

    else{
        c <- c + freq[i]*
          (log( pnorm(theta[2]*TL[i+1]-theta[1])-pnorm(theta[2]*TL[i]-theta[1])))
    }#end else
  }#end for

  L <- -(a+b+c)
  
  return(L)
}
################################################################################


### exact_mle FUNCTION -------------------------------------------------------------

#' Maximizing the Likelihood Estimation.
#'
#' @param TL Lower end of the intervals.
#' @param freq Frequency on each interval.
#' @param intial_mu Initial numeric value for mu (mean).
#' @param intial_sigma Initial numeric value for sigma.
#' @return Return the result of the log likelihood for the mean (mu) and sigma 
#'  (variance) estimates.
#' @export
#' @examples
#'  simdataaaa = univ_Simul(ncol_matrix=1,
#'                          n=50,
#'                          nclass = 10,
#'                          mean = 68,
#'                          sd = 1.80,
#'                          fr_breaks=c(62,64,66,68,70,72,74,76,78))
#'  
#'  
#'  TL<- simdataaaa$simul_data[,1,1]
#'  freq<- simdataaaa$simul_data[,3,1]
#'  
#'  estimates_simul = exact_mle(TL=TL,
#'                              freq=freq,
#'                              intial_mu = (67/2),
#'                              intial_sigma = (1/2))
#'  
#'  estimates_simul

exact_mle <- function(TL = c(-Inf,62,64,66,68,70,72,74,76,78),
                      freq = c(0,1,9,16,15,8,1,0,0,0),
                      intial_mu = (67/2),
                      intial_sigma = (1/2)){
  
  res <- stats::optim(c(intial_mu,intial_sigma),
                      fn=Logll,
                      TL=TL,
                      freq=freq,
                      method="L-BFGS-B",
                      lower=c(0.02,0.2),
                      upper=c(180,2))
  
  mu_estimate<- base::round(res$par[1]/res$par[2],4)
  sigma_estimate<- base::round(1/res$par[2],4)
  sigma_squared_estimate = sigma_estimate^2
  
  final_list <- base::list("mu_estimate" = mu_estimate,
                           "sigma_estimate" = sigma_squared_estimate)
  
  return(final_list)
  
}
################################################################################









