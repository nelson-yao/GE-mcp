library(mvtnorm)

#n = 100
#p = 1000
#q = 5
 
generate <- function(k,mag2,n,p,q,quan)
{
 
  covar.z <- matrix(rep(0,p*p),nrow=p)

  for ( i in 1:p)
   for ( j in 1:p)
   {
    if ( abs(i-j) == 1)
     covar.z[i,j]=0.33
   }

  for (i in 1:p)
  {
    covar.z[i,i]=1
  }

# system.time(z <- rmvnorm(n = n, mean = rep(0, p), covar.z))
# system.time(x <- rmultnorm(n,rep(0,p),covar) )
  Z <- rmvnorm(n = n, mean = rep(0, p), covar.z)
  
  covar.x <- matrix(rep(0,(q-1)*(q-1)),nrow=(q-1))

  for ( i in 1:(q-1))
    for ( j in 1:(q-1))
    covar.x[i,j]=mag2^(abs(i-j))

  X <- matrix(rep(0,n*q),nrow=n)
  #X <- rmvnorm(n = n, mean = rep(0, q), covar.x)
  X[,1:(q-1)] <- rmvnorm(n = n, mean = rep(0, (q-1)), covar.x)
  X[,q] <- rbinom(n, 1, 0.6)
  
  index=seq(1,p)
  fun <- function(i,Z,X)
  {
    V.i=Z[,i]*X
    list(V.i=V.i)
  }
  tmp=lapply(index,fun,Z,X)
  V=matrix(unlist(tmp),nrow=n) 

  W <- matrix(rep(0,n*p*(q+1)), n);
  for ( j in 1:p)
  {
    W[,(j-1)*(q+1)+1] <- Z[,j];
    W[,((j-1)*(q+1)+2):(j*(q+1))] <- V[,((j-1)*q+1):(j*q)];
  }

  I <- diag(q);
  b <- numeric(p*(q+1))
  for ( j in 1:(2*q))
  {
     tmp <- j%%q;
     if ( tmp == 0) tmp <- q;
     b[((j-1)*(q+1)+1): (j*(q+1))]=c(1,rep(1,q)*I[,tmp]);
     #b[((j-1)*(q+1)+1): (j*(q+1))]=c(runif(1,0.5,1.5),runif(5,0.5,1.5)*I[,tmp]);
     #b[((j-1)*(q+1)+1): (j*(q+1))]=c(1,2*I[,tmp]);
  }

  #alpha <- runif(5,0.5,1.5)
  alpha <- rep(1,q)
  res <-rnorm(n = n, 0 , k) 
  T1 <- X%*%alpha+W%*%b+res
  #T1 <- W%*%b+res
  
  C <- runif(n,min=0,max=quantile(exp(T1),c(quan)))
  Y <- pmin(exp(T1), C)
  D <- as.numeric(exp(T1) <= C)
  Xs=X[order(Y),]
  Ws=W[order(Y),]
  D=D[order(Y)]
  Y=sort(Y)

  weight <- numeric(n)
  weight[1]=D[1]/n
  for ( i in 2:n )
  {
    tmp = 1
    for ( j in 1: (i-1) )
      tmp = tmp*((n-j)/(n-j+1))^D[j]
     
    weight[i]=D[i]/(n-i+1)*tmp
  }

  #snr <- (t(x%*%beta)%*%(x%*%beta)/n/k^2)^0.5
  list(X=X,Z=Z,V=V,Y=Y,W=W, Xs=Xs, Ws=Ws, D=D, C=C, weight=weight, T1=T1, b=b)
}

#source("/space/jliu20/GE_interaction/generate.surv.band1.r")
