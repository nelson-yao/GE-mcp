dyn.load("/space/jliu20/GE_interaction/sparse_GE.so")

GE_interact <- function(Y,X,W,lambda1,lambda2,gamma1,gamma2)
{
  #Z is p genomic risk factor; X is q clinical/Environmental factor; V is pq interaction terms.
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(W)/(q+1)
  
  group <- rep(q+1,p)
  totalp <- sum(group)
  
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
        
  old.beta <- numeric(p*(q+1));

  W.t <- W
  
  diff <- 1 
  epsilon <- 1E-8
  param <- c(n, p, totalp,200) 
  iter <- numeric(1)
  alpha <- numeric(q)
  repeat 
  {
    Y.p <- Y - X%*%alpha;
    fit <-  .C("sparse_GE", y=as.double(Y.p),x=as.double(t(W.t)),G=as.integer(group),param=as.integer(param),
                          tuning=as.double(c(lambda1,lambda2,gamma1,gamma2)),
                          epsilon=as.double(epsilon),beta=as.double(old.beta),diff=as.double(diff))
    
    r_w <- Y - W.t%*%fit$beta;
    alpha <- solve(t(X)%*%X)%*%t(X)%*%r_w;
    diff <- max(abs(old.beta-fit$beta));
    old.beta  <- fit$beta;     
    iter <- iter + 1;              
    if ( diff < 1E-5 |  iter >= 50 )
      break;
    
  }


  list(alpha=alpha,beta=fit$beta,diff=diff,W.t=W.t)
}


#source("/space/jliu20/GE_interaction/CallC_iterate.r")