#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
arma::cube gibbsMessfull(int iter1, arma::mat Ainv, arma::mat expo, double kappa, arma::mat X, int kdim, double theta0, int nn, int iter)
{
	/*input: iter1=burn-in of Gibbs sampler, Ainv: sigma^2*Ainv is the covariance
	 matrix of the multivariate Normal distribution of beta, expo=matrix of vectors expm(rho*D) for all rho, 
	 * where D is the spatial weights matrix, kappa=shape parameter of inverse gamma distribution of sigma^2
	  X=predictor variables, kdim=number of predictor variables, theta0=scale parameter of the prior distribution of sigma^2, 
	   nn=number of samples of the outcome variable, iter=number of draws from distribution of rho
	   output: matrix containing samples from the joint posterior distribution of beta and sigma^2*/
	   arma::cube MM(499,kdim+1,iter);
	   arma::mat exponential(nn,1);
	   arma::mat helpVec, factor,thetaMat,betaMeanTemp;
	   arma::mat M(iter1,kdim+1);
	   int i,j,k,l;
	   double a;
	   for ( k=0;k<iter-1;k++)
	   {
		   
	for (l=0;l<nn;l++)	
	{
		exponential(l,0)=expo(l,k);
	}
betaMeanTemp=Ainv*(X.t()*exponential);
	for (i=0;i<kdim;i++)
	{
		M(0,i)=betaMeanTemp(i);
	}
	M(0,kdim)=0.5;
	for (i=0;i<(iter1-1);i++)
	{
	helpVec=betaMeanTemp.t()+arma::randn(1,kdim) * arma::chol(Ainv*M(i,kdim));
	for (j=0;j<kdim;j++)
	{
	M(i+1,j)=helpVec(0,j);}
	factor=exponential-X*helpVec.t();
	/*Rcpp::Rcout<<theta;*/
	thetaMat=factor.t()*factor;
	double theta=(0.5*as_scalar(thetaMat)+theta0);
	a=rgamma(1,kappa,1/theta)(0);
	M(i+1,kdim)= 1/a;
	}
	MM.slice(k)=M.submat(iter1-500,0,iter1-2,kdim);
}
	return MM;
}
	