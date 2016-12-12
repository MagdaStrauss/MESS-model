#This R file includes a function for Bayesian estimation of the MESS error model (with splines).
#It also contains code for an application to house price data.


require(R.matlab)
require(matrixcalc)
require(expm)
require(MHadaptive)
require(graphics)
require(coda)
require(MCMCpack)
require(spdep)
require(CARBayes)
require(splines)
require(Rcpp)
require(RcppArmadillo)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("MESSGibbsExtMod.cpp")

ExtendModMargRho<-function(X,y,D,iter)
  #This function samples from the marginal posterior distribution of the spatial
  #parameter rho in the MESS error model (with or without splines). 
  
  #If splines are used, then they are included in X.
  
  #input: X...predictors, y...outcome, D...weights matrix, 
  #iter=number of samples to be drawn
  #output: samples from the marginal posterior distribution of rho
{
  kappa0=0.001
  theta0=0.001 
  kdim<-dim(X)[2]
  Cinv=10^(-4)*diag(kdim)
  nn=length(y)
  kappa<-kappa0+nn/2
  logPdf<-function(rho)
  {
    expon<-expm(rho*D)
    M<-expon%*%X
    H<-diag(nn)-M%*%chol2inv(chol(t(M)%*%M))%*%t(M)
    bhat<-chol2inv(chol(t(M)%*%M))%*%t(M)%*%expon%*%y
    A=t(M)%*%M+Cinv
    Ainv=chol2inv(chol(A))
    btilde<-Ainv%*%(t(M)%*%M%*%bhat)
    Q1<-t(btilde)%*%Cinv%*%btilde
    M<-expon%*%X
    Q2<-t(bhat-btilde)%*%(t(M)%*%M)%*%(bhat-btilde)
    factor<-H%*%expon%*%y
    PP<-t(factor)%*%factor
    Htilde=
      output<-(-0.5*(nn+2*kappa0))*log(2*theta0+PP+Q1+Q2)+log(det(Ainv))
    return(output)
  }
  E<-function(rho)
  {
    expon<-expm(rho*D)
    M<-expon%*%X
    H<-diag(nn)-M%*%chol2inv(chol(t(M)%*%M))%*%t(M)
    factor<-H%*%expon%*%y
    PP<-t(factor)%*%factor
    return(PP)
  }
  
  x<-optimize(E, interval=c(-5, 1), maximum=FALSE)
  samples<-Metro_Hastings(li_func=logPdf, pars=x$minimum,  prop_sigma = NULL,
                          par_names = NULL, iterations = 5000+10*iter, burn_in = 5000,
                          adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE)
  samples<-mcmc_thin(samples, thin = 10) 
  samples<-mcmc(samples$trace)
  densplot(samples)
  autocorr.plot(samples)
  return(samples)}
##############
ExtendModGibbs<-function(X,y,D,iter,iter1,samples)
{
  #This function samples from the marginal posterior distribution of the spatial parameter rho.
  
  #input: X=independent variables, y=outcome, D=weights matrix, iter=number of samples of rho used in the Gibbs sampler 
  #iter1=burn-in of Gibbs sampler, samples=samples from the marginal posterior of rho
  #output: samples from the joint posterior of beta and sigma^2, DIC
  
  kappa0=0.001
  theta0=0.001 
  kdim<-dim(X)[2]
  samples<-samples[1:iter]
  Cinv=10^(-4)*diag(kdim)
  nn=length(y)
  kappa<-kappa0+nn/2
  MESSExt<-c()
  MESSExt$beta<-array(data=NA,dim=c(kdim,iter))
  MESSExt$sigma2<-matrix(data=NA,nrow=1,ncol=iter)
  for(i in seq(1,iter,by=1))
  {
    print(i)
    expo<-expm(samples[i]*D)
    expoX<-expo%*%X
    expoY<-expo%*%y
    covas<-chol2inv(chol(t(expoX)%*%expoX+Cinv))
    betaMeanTemp=covas%*%t(expoX)%*%expoY
    AA<-gibbsMessExtMod(iter1+1,betaMeanTemp,covas,expo,kappa,X,y,kdim,theta0)
    MESSExt$beta[, i]<-AA[(iter1+1),1:kdim]
    MESSExt$sigma2[i]<-AA[(iter1+1),kdim+1]}
  MESSExt$rho<-samples[1:(iter-1)]
  #### DIC
  MESSLikeE<-function(pars)
  {
    lp=length(pars)
    beta<-matrix(data=pars[2:(lp-1)],nrow=lp-2,ncol=1)
    A<-expAtv(D,y-X%*%beta,pars[1])$eAtv
    output<-(-2)*(-nn/2*log(2*pi)-nn/2*log(pars[lp])-1/(2*pars[lp])*t(A)%*%A)
    return(output)
  }
  vec1<-matrix(data=MESSExt$rho,nrow=1,ncol=iter-1)
  helpMat<-rbind(vec1,MESSExt$beta[,1:(iter-1)],matrix(data=MESSExt$sigma2[(1:iter-1)],nrow=1,ncol=iter-1))
  helpMatrix<-apply(helpMat,2,MESSLikeE)
  DIC.av<-mean(helpMatrix)
  rho.mean<-mean(vec1)
  beta.mean<-rowMeans(MESSExt$beta[,1:(iter-1)])
  dim(beta.mean)<-c(1,kdim)
  sigma2.mean<-mean(MESSExt$sigma2[(1:iter-1)])
  MESSExt$pD<-DIC.av-MESSLikeE(cbind(rho.mean,beta.mean,sigma2.mean))
  MESSExt$DIC<-MESSExt$pD+DIC.av
  MESSExt$pV<-0.5*var(helpMatrix)
  MESSExt$DICV<-MESSExt$pV+DIC.av
  return(MESSExt)
}


#####House price example including splines
load("spatialhousedata.rda")
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
M.nb <- poly2nb(spatialhousedata, row.names = rownames(spatialhousedata@data))
M.list <- nb2listw(M.nb, style = "B")
M.mat <- nb2mat(M.nb, style = "W")

#exporting the housing data to Matlab
writeMat("houseX.mat",A=house.X)
writeMat("houseMatrix.mat",W=M.mat)
writeMat("log.house.price.mat",B=log.house.price)
#saving as R data files
save(house.X,file="house.X")
save(log.house.price,file="log.house.price")
save(M.mat,file="houseMatrix")

cc<-coordinates(spatialhousedata)
spline1<-ns(cc[,1],df=3)
spline2<-ns(cc[,2],df=3)
X.extra<-matrix(data=NA,nrow=270,ncol=15)
X.extra[,1:3]<-spline1
X.extra[,4:6]<-spline2
k=7
for (i in seq(1,3,1))
{
  for(j in seq(1,3,1))
  {
    X.extra[,k]=spline1[,i]*spline2[,j]
    k=k+1
  }
}
house.XX<-cbind(house.X,X.extra)
SplinesRho1<-ExtendModMargRho(house.XX,log.house.price,M.mat,5000)
SplinesVersion1<-ExtendModGibbs(house.XX,log.house.price,M.mat,5000,4000,SplinesRho1)
### only with individual splines, excluding the tensor product
house.X2<-house.XX[,1:14]
spline1<-ns(cc[,1],df=5)
spline2<-ns(cc[,2],df=5)
X.extra<-matrix(data=NA,nrow=270,ncol=10)
X.extra[,1:5]<-spline1
X.extra[,6:10]<-spline2
house.X3=cbind(house.X,X.extra)
SplinesRho3<-ExtendModMargRho(house.X3,log.house.price,M.mat,5000)
SplinesVersion3<-ExtendModGibbs(house.X3,log.house.price,M.mat,5000,4000,SplinesRho3)

