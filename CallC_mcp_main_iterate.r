dyn.load("/space/jliu20/GE_interaction/gmcp.so")

mcp <- function(X,Z,Y,lambda,gamma)
{
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(Z)
  
  group <- rep(1,p)
  totalp <- sum(group)
  
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
  
  fun <- function(i,Z)
  {
    sigma = t(Z[,pos_s[i]:pos_e[i]])%*%Z[,pos_s[i]:pos_e[i]]/n
    R = chol(sigma)  #cholesky decomp always exists for X'X/n where dj < n
    Z.t = Z[,pos_s[i]:pos_e[i]]%*%solve(R)
    return(list(R=R,Z.t=Z.t))
  }
  result <- lapply(seq(1,totalp),fun,Z) 

  fun <- function(j,result)
  {
    result[[j]]$Z.t
  }
  Z.t <- matrix(unlist(lapply(seq(1,totalp),fun,result)),nrow=n)

  old.beta <- numeric(p);

  diff <- 1 
  epsilon <- 1E-8
  param <- c(n, totalp, totalp,200,1)  
  iter <- numeric(1)
  alpha <- solve(t(X)%*%X)%*%t(X)%*%Y
  repeat 
  {
    Y.p <- Y - X%*%alpha;
    fit <- .C("GMCP", y=as.double(Y.p),x=as.double(t(Z.t)),G=as.integer(group),param=as.integer(param),
                      lambda=as.double(lambda),gamma=as.double(gamma),epsilon=as.double(epsilon),
                      beta=as.double(old.beta),diff=as.double(0))

    #fit <- .C("sparse_GE", y=as.double(Y.p),x=as.double(t(Z.t)),G=as.integer(group),param=as.integer(param),
    #                      tuning=as.double(c(lambda1,lambda2,gamma1,gamma2)),
    #                      epsilon=as.double(epsilon),beta=as.double(old.beta),diff=as.double(diff))
    
    r_w <- Y - Z.t%*%fit$beta;
    alpha <- solve(t(X)%*%X)%*%t(X)%*%r_w;
    diff <- max(abs(old.beta-fit$beta));
    old.beta  <- fit$beta;     
    iter <- iter + 1;              
    if ( diff < 1E-5 |  iter >= 50 )
      break;
    
  }

  beta <- numeric()
  b = fit$beta
  for ( j in 1:totalp )
  {
     beta <- c(beta,solve(result[[j]]$R,b[pos_s[j]:pos_e[j]]))
  }

  list(alpha=alpha,beta=beta,b=b,diff=fit$diff,result=result,Z.t=Z.t)
}

#source("/space/jliu20/GE_interaction/CallC_mcp_main_iterate.r")
