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

make.X0 <- function(n, rho, p0){
  X0 <- (mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p0))), sigma = make.SIGMA(rho, p0)))
}


make.X1 <- function(n, rho, p0, p1){
  X0 <- make.X0(n, rho, p0)
  X1 <- mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p1))), sigma = make.SIGMA(rho, p1))
  
  ((diag(1,n) - X0 %*% solve(t(X0) %*% X0) %*% t(X0)) %*%  X1)
  
}

make.B0 <- function(n, p0, R2_0, sigma2, X0){
  
  Jn <- diag(n) - (matrix(rep(1,n),n,1) %*% t(matrix(rep(1,n),n,1)))/n
  K0 <- t(matrix(rep(1,p0),p0,1)) %*% t(X0) %*% Jn %*%  X0 %*% matrix(rep(1,p0),p0,1)
  num <- sigma2 *((n-1)*R2_0 - (p0 - 1))
  den <- K0 *(1-R2_0)
  
  return(sqrt(num/den))
  
}  

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

make.y.under.H0 <- function(n, p0, sigma2,R2_0, X0){
  
  (X0 %*%  (matrix(rep(1,p0), p0 ,1) %*% make.B0(n, p0, R2_0, sigma2, X0)) + rnorm(n,0, sigma2))
  
}

make.y.under.H1 <- function(n, p0, p1, sigma2,R2_0, R2_1, X0,X1){
  
  (X0 %*%  (matrix(rep(1,p0), p0 ,1) %*% make.B0(n, p0, R2_0, sigma2, X0))
   + X1 %*%  (matrix(rep(1,p1), p1 ,1) %*% make.B1(n, p0, p1, R2_0,R2_1, sigma2, X0,X1))
   + rnorm(n,0, sigma2))
  
}

cutoff_finder_normal = function(p1, n, delta)
{
  f1 = function(T) pchisq(T, df=p1, lower.tail = FALSE) - pnorm( (T-p1-n*delta)/sqrt(2*(p1+2*n*delta)) )
  uniroot(f1, lower=0, upper=1e8, extendInt="yes", maxiter = 1e9, tol =  .Machine$double.eps^0.75)$root
  
}



