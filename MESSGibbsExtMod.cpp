#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
arma::mat gibbsMessExtMod(int iter1, arma::colvec betaMeanTemp, arma::mat Ainv, arma::mat exponential, double kappa, arma::mat X,arma::colvec y,int kdim, double theta0)
{
	/*input: iter1=burn-in of Gibbs sampler, betaMeanTemp=mean of multivariate Normal distribution of beta, Ainv: sigma^2*Ainv is the covariance
	 matrix of the multivariate Normal distribution of beta, exponential=exp(rho*D)*y, kappa=shape parameter of inverse gamma distribution of sigma^2
	  X=predictor variables, y=outcome, kdim=number of predictor variables, theta0=scale parameter of the prior distribution of sigma^2, 
	   output: matrix containing samples from the joint posterior distribution of beta and sigma^2*/
	arma::rowvec betaMeanTempt=betaMeanTemp.t();
	arma::mat M(iter1,kdim+1);
	int i;
	int j;
	double theta,a;
	arma::mat factor, thetaMat, helpVec;
	for (i=0;i<kdim;i++)
	{
		M(0,i)=betaMeanTemp(i);
	}
	M(0,kdim)=0.5;
	for (i=0;i<(iter1-1);i++)
	{
	helpVec=betaMeanTempt+arma::randn(1,kdim) * arma::chol(Ainv*M(i,kdim));
	for (j=0;j<kdim;j++)
	{
	M(i+1,j)=helpVec(0,j);}
	factor=exponential*(y-X*helpVec.t());
	/*Rcpp::Rcout<<theta;*/
	thetaMat=factor.t()*factor;
	theta=(0.5*as_scalar(thetaMat)+theta0);
	a=Rcpp::rgamma(1,kappa,1/theta)(0);
	M(i+1,kdim)= 1/a;
	}
	return M;
}