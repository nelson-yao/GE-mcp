home="~/Dropbox/workfolder/ABstudy/Codes-Liu/"
dyn.load(paste0(home,"gmcp.so"))


mcp <- function(X,Z,Y,lambda,gamma)
{
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(Z)
  
  group <- rep(1,p)
  totalp <- sum(group)
  
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
  
  Z.t <- Z

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

  list(alpha=alpha,beta=fit$beta,diff=fit$diff,Z.t=Z.t)
}

#source("/space/jliu20/GE_interaction/CallC_mcp_iterate_wo.r")
