# Creates the Variance-covariance matrix

make.SIGMA <- function(rho,dimension){
  SIGMA = matrix(NA,dimension,dimension)
  for(i in 1:dimension){
    for(j in 1:dimension){
      if(i != j){
        a <- sign(rnorm(1,0,1))*rnorm(1,0.1, 0.01)
        if(a>0) {SIGMA[i,j] = rho + a}
        else if (a <0)  {SIGMA[i,j] = rho}
      }
      else if (i == j){
        SIGMA[i,j] = 1
      }
    }
  }
  (SIGMA + t(SIGMA))/2
}

# Design matrix for baseline covariates and confounders

make.X0 <- function(n, rho, p0){
  X0 <- (mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p0))), sigma = make.SIGMA(rho, p0)))
}

# Design matrix for exposures or exposure mixture 

make.X1 <- function(n, rho, p0, p1){
  X0 <- make.X0(n, rho, p0)
  X1 <- mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p1))), sigma = make.SIGMA(rho, p1))
  
  ((diag(1,n) - X0 %*% solve(t(X0) %*% X0) %*% t(X0)) %*%  X1)
}

# Beta coefficents for baseline covariates such that the $R^2$ is fixed

make.B0 <- function(n, p0, R2_0, sigma2, X0){
  
  Jn <- diag(n) - (matrix(rep(1,n),n,1) %*% t(matrix(rep(1,n),n,1)))/n
  K0 <- t(matrix(rep(1,p0),p0,1)) %*% t(X0) %*% Jn %*%  X0 %*% matrix(rep(1,p0),p0,1)
  num <- sigma2 *((n-1)*R2_0 - (p0 - 1))
  den <- K0 *(1-R2_0)
  
  return(sqrt(num/den))
  
}  

# Beta coefficents for Exposures such that the $R^2$ is fixed for the larger model given the beta coefficients for the smaller model


make.B1 <- function(n, p0, p1, R2_0, R2_1, sigma2,X0,X1){
  
  Jn <- diag(n) - (matrix(rep(1,n),n,1) %*% t(matrix(rep(1,n),n,1)))/n
  B0 <- make.B0(n, p0, R2_0, sigma2, X0)
  
  K00 <- t(as.matrix(rep(1,p0)))%*% t(X0) %*% Jn %*% X0 %*% as.matrix(rep(1,p0))
  K11 <- t(as.matrix(rep(1,p1)))%*% t(X1) %*% Jn %*% X1 %*% as.matrix(rep(1,p1))
  K10 <- t(as.matrix(rep(1,p1)))%*% t(X1) %*% Jn %*% X0 %*% as.matrix(rep(1,p0))
  K01 <- t(as.matrix(rep(1,p0)))%*% t(X0) %*% Jn %*% X1 %*% as.matrix(rep(1,p1))
  
  A <- K11 *(1 - R2_1)
  B <-  B0 * (K10 + K01) *( 1 - R2_1)
  C <- (B0^2 * K00 + sigma2 * (p0 + p1 - 1) - R2_1 *((n-1)*sigma2 + B0^2 * K00))  
  return((- B + sqrt(B^2 - 4*A*C))/(2*A))
  
}

# Generate y under the null

make.y.under.H0 <- function(n, p0, sigma2,R2_0, X0){
  
  (X0 %*%  (matrix(rep(1,p0), p0 ,1) %*% make.B0(n, p0, R2_0, sigma2, X0)) + rnorm(n,0, sigma2))
}

# Generate y under the alternative

make.y.under.H1 <- function(n, p0, p1, sigma2,R2_0, R2_1, X0,X1){
  
  (X0 %*%  (matrix(rep(1,p0), p0 ,1) %*% make.B0(n, p0, R2_0, sigma2, X0))
   + X1 %*%  (matrix(rep(1,p1), p1 ,1) %*% make.B1(n, p0, p1, R2_0,R2_1, sigma2, X0,X1))
   + rnorm(n,0, sigma2))
  
}

# Numerically estimate T, the calibrated cutoff with specified $p_1$, n and $\delta$

cutoff_finder_normal = function(p1, n, delta)
{
  f1 = function(T) pchisq(T, df=p1, lower.tail = FALSE) - pnorm( (T-p1-n*delta)/sqrt(2*(p1+2*n*delta)) )
  uniroot(f1, lower=0, upper=1e8, extendInt="yes", maxiter = 1e9, tol =  .Machine$double.eps^0.75)$root
  
}

# Numerical estimate of rate of change of T as n \to \infty

numerical_rate_normal = function(p1, n,  delta){
  
  n1 = 3e4; n2 = 3e4+1   
  (cutoff_finder_normal(p1, n2, delta) - cutoff_finder_normal(p1, n1, delta))/(n2 - n1)
  
}

# Calculates type 2 error  

beta_func <- function(p1, n, delta, epsilon){
  T = cutoff_finder_normal(p1, n, delta)
  pnorm((T - p1 - n*epsilon)/sqrt(2*(p1 + 2*n*epsilon)))
}




