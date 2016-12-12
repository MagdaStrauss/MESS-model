#This file includes the functions used for our Bayesian inference
#of the MESS lag model. It also contains a function to compute the CPO
#and an application to house price data.

require(R.matlab)
require(matrixcalc)
require(expm)
require(MHadaptive)
require(graphics)
require(coda)
require(MCMCpack)
require(sp)
require(spdep)
require(Rcpp)
require(RcppArmadillo)
require(RANN)
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("MESSGibbsfull.cpp")
#####################


cpo<-function(rho,beta,sigma,D,y,X)
{
	#Computes the CPO for the MESS lag model
	
	#input: rho, beta, sigma ... parameters of the MESS lag model
	# D.... spatial weights matrix
	#y...outcome
	#X...predictors
	
	
  A<-expAtv(D,y-X%*%beta,rho)$eAtv
  output=0*y
  for (j in 1:length(output)){
    output[j]<-(-1/2*log(2*pi)-1/2*log(sigma)-1/(2*sigma)*t(A[j])%*%A[j])}
  #print(output)
  return(output)
}


RhoMargPost<-function(X,y,D,iter)
{
	
#This function samples from the marginal posterior
#of the spatial parameter rho for the MESS lag model.
#It uses the function 'expAtv' from the 'expm' package
#to compute the action of a matrix exponential on a vector.

#input: #X...predictors
	#y...outcome
	#D...spatial weights matrix
	#iter...number of draws from distribution of rho
#output: samples from marginal posterior of rho

  kappa0=0.001
  theta0=0.001 
  kdim<-dim(X)[2]
  C=10^4*diag(kdim)
  C.inv<-diag(kdim)*10^(-4)
  A=t(X)%*%X+C.inv
  A.inv<-chol2inv(chol(A))
  nn=length(y)
  kappa<-kappa0+nn/2
  H<-diag(nn)-X%*%chol2inv(chol(t(X)%*%X))%*%t(X)
  logPdf<-function(rho)#marginal log-density of rho
  {
    XX<-t(X)%*%X
    XX.inv<-chol2inv(chol(XX))
    expon=expAtv(D,y,rho)$eAtv
    bhat<-XX.inv%*%t(X)%*%expon
    btilde<-A.inv%*%(t(X)%*%X%*%bhat)
    Q1<-t(btilde)%*%C.inv%*%(btilde)
    Q2<-t(bhat-btilde)%*%(t(X)%*%X)%*%(bhat-btilde)
     factor<-H%*%expAtv(D,y,rho)$eAtv
     PP<-t(factor)%*%factor
    output<-(-0.5*(nn+2*kappa0))*log(2*theta0+PP+Q1+Q2)
    return(output)
  }
  E<-function(rho)
  {expon=expAtv(D,y,rho)$eAtv
   ouput<-matrix(data=NaN,1,1)
   factor<-H%*%expon
   output<-t(factor)%*%factor
   return(output)}
  x<-optimize(E, interval=c(-5, 1), maximum=FALSE)
  samples<-Metro_Hastings(li_func=logPdf, pars=x$minimum,  prop_sigma = NULL,
                          par_names = NULL, iterations = 10*iter+5000, burn_in = 5000,
                          adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE)
  samples<-mcmc_thin(samples, thin = 10) #remove autocorrelation
  samples<-mcmc(samples$trace)
  densplot(samples)
  autocorr.plot(samples)
 return(samples)}
 

#####################
#Gibbs sampler using C++ function
  MESSAuto2<-function(X,y,D,iter,iter1,samples)
  {
  	#This function samples from the posteriors 
  	#of the non-spatial MESS paramters
  	#given the samples from the posterior
  	#of rho
  	#iter1=burn-in of Gibbs sampler
  	#iter = number of samples from the marginal of rho
  	#samples= samples from the marginal posterior of rho
  kappa0=0.001
  theta0=0.001 
  kdim<-dim(X)[2]
  C=10^4*diag(kdim)
  C.inv<-diag(kdim)*10^(-4)
  A=t(X)%*%X+C.inv
  A.inv<-chol2inv(chol(A))
  C.inv<-diag(kdim)*1/10000
  nn=length(y)
  print(nn)
  kappa<-kappa0+nn/2
  expofun<-function(rho,D,y)
  {return(expm::expAtv(D,y,rho)$eAtv)}
  samples=samples[-iter]
  samples=matrix(data=samples[1:iter-1],nrow=iter-1,ncol=1)
  expo=apply(samples,1,expofun,D=D,y=y)
   AA<-gibbsMessfull(iter1+499,A.inv,expo,kappa,X,kdim, theta0,nn,iter)
  MESSAuto<-c()
  MESSAuto$rho<-samples[1:(iter-1)]
  MESSAuto$beta<-AA[,1:kdim,1:(iter-1)]
  MESSAuto$sigma2<-AA[,kdim+1,1:(iter-1)]
  
  MESSLike<-function(pars)
{
  #input: pars=rho, beta, sigma^2....parameters of the MESS model
  #output: (-2)*logarithm of MESS-Likelihood (required to compute DIC)
  
  lp<-length(pars)
  beta<-matrix(data=pars[2:(lp-1)],nrow=lp-2,ncol=1)
  A<-expAtv(D,y,pars[1])$eAtv-X%*%beta
  output<-(-2)*(-nn/2*log(2*pi)-nn/2*log(pars[lp])-1/(2*pars[lp])*t(A)%*%A)
  return(output)
}
 
vec1<-matrix(data=MESSAuto$rho,nrow=1,ncol=iter-1)
  helpMat<-rbind(vec1,AA[499,,1:(iter-1)])
  helpMatrix<-apply(helpMat,2,MESSLike)
  DIC.Av<-mean(helpMatrix)
  rho.mean<-mean(MESSAuto$rho)
  beta.mean<-rowMeans(MESSAuto$beta[499,,])
  sigma2.mean<-mean(MESSAuto$sigma2[499,])
  dim(beta.mean)<-c(length(beta.mean),1)
  pD<-DIC.Av-MESSLike(pars=c(rho.mean,beta.mean,sigma2.mean))
  DICAuto<-pD+DIC.Av
  pV<-var(helpMatrix)/2
  DICAutoV<-pV+DIC.Av
  MESSAuto$DIC<-DICAuto
  MESSAuto$pD<-pD
  MESSAuto$DICV<-DICAutoV
  MESSAuto$pV<-pV
  return(MESSAuto)
}

#house price application

require(CARBayes)
load("spatialhousedata.rda")#data provided by the CARBayes package
#(Scottish Neighbourhood Statistics)
#with row-stochastic contiguity matrix
M.nb <- poly2nb(spatialhousedata, row.names = rownames(spatialhousedata@data))
M.list <- nb2listw(M.nb, style = "B")
M.mat <- nb2mat(M.nb, style = "W")
log.house.price<-log(spatialhousedata@data$price)
house.rooms<-spatialhousedata@data$rooms
house.crime.log<-log(spatialhousedata@data$crime)
house.sales<-spatialhousedata@data$sales
house.driveshop.log<-log(spatialhousedata@data$driveshop)
house.type<-spatialhousedata@data$type
house.type.semi<-as.numeric(house.type=="semi")
house.type.flat<-as.numeric(house.type=="flat")
house.type.terrace<-as.numeric(house.type=="terrace")
intercept<-matrix(data=1,nrow=length(log.house.price), ncol=1)
house.X<-c(intercept,house.rooms,house.crime.log,house.sales,house.driveshop.log,house.type.semi,house.type.flat,house.type.terrace)
dim(house.X)<-c(270,8)
dim(log.house.price)<-c(270,1)
#exporting the housing data to Matlab
writeMat("houseX.mat",A=house.X)
writeMat("houseMatrix.mat",W=M.mat)
writeMat("log.house.price.mat",B=log.house.price)
#saving as R data files
save(house.X,file="house.X")
save(log.house.price,file="log.house.price")
save(M.mat,file="houseMatrix")
RhoMarg<-RhoMargPost(house.X,log.house.price,M.mat,5000)
SamplesAll<-MESSAuto2(X=house.X,y=log.house.price,D=M.mat,iter=5000,iter1=4000,samples=RhoMarg)
#save(SamplesAll, file="MESSrowstochCont")
####

#with 7-nearest neighbours matrix
coords<-coordinates(spatialhousedata)
M7<-knearneigh(coords,k=7)
M7<-knn2nb(M7)
M7<-nb2mat(M7,style="W")
RhoMarg7<-RhoMargPost(house.X,log.house.price,M7,5000)
SamplesAll7<-MESSAuto2(house.X,log.house.price,M7,5000,4000,RhoMarg7)
#with binary contiguity matrix
Mbin.mat <- nb2mat(M.nb, style = "B")
RhoBin<-RhoMargPost(house.X,log.house.price,Mbin.mat,5000)
SampleAllBin<-MESSAuto2(house.X,log.house.price,Mbin.mat,5000,4000,RhoBin)



