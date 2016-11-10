#.libPaths("~/rpacks")
library(snow)

lammax <- function(x,y,group)
{
  p <- length(group)
  n <- dim(x)[1]
  totalp <- ncol(x)

  #calculate lambda max
  end = cumsum(group)
  start = end -group + 1
  inner <- rep(0,p)
  for ( j in 1:p )
  {
    inner[j]=sum(x[,start[j]:end[j]]*as.vector(y)/n) 
  }
  lambda.max = max(inner)
  list(lambda.max = lambda.max)
}

sse.evaluate <- function(i,X,Y,W,fold.d,group,lambda,gamma)
{
  delete  <- unlist(fold.d[i])
  X.train <- X[-delete,]
  Y.train <- Y[-delete]
  W.train <- W[-delete,]
  
  X.valid <- X[delete,]
  Y.valid <- Y[delete]
  W.valid <- W[delete,]

  source("/space/jliu20/GE_interaction/CallC_mcp_iterate.r")
  #source("/space/jliu20/GE_interaction/CallC_mcp_iterate.wo.r")
  fit <- mcp(X.train,W.train,Y.train,lambda,gamma)
  #fit <- GE_interact(Y.train,X.train,W.train,lambda1,lambda2,gamma,gamma)
  
  sse <- sum((Y.valid - X.valid%*%fit$alpha - W.valid%*%fit$beta)^2)
  #sse <- sum((y.valid -fit$beta0- x.valid%*%fit$beta)^2)
  list(sse=sse)
}

parCv <- function(cl,lambda.seq,X,Y,W,group,fold.d,n.fold,gamma)
{
  step <- length(lambda.seq)
  PE <- numeric(step)
  index <- seq(1:n.fold)

  for( i in 1:step)
  {
    lambda=lambda.seq[i]
    #v <- lapply(index, sse.evaluate, X,Y,W,fold.d,group,lambda,gamma)
    v <- clusterApply(cl, index, sse.evaluate, X,Y,W,fold.d,group,lambda,gamma)
    PE[i] <- sum(unlist(v))/n.fold
  }
  list(PE=PE)
}

#cv.optim <- function(X,Y,W,group,n.fold,epsilon1,epsilon2,alpha,n.step,gamma,cl)
cv.optim <- function(X,Y,W,n.fold,epsilon,n.step,gamma,cl)
{
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(W)/(q+1)
  
  group <- rep(1,(q+1)*p)
  totalp <- sum(group)
  
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1

  fold.d = splitList(sample(c(1:n),n),n.fold)    

  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1

  fun <- function(i,W)
  {
    sigma = t(W[,pos_s[i]:pos_e[i]])%*%W[,pos_s[i]:pos_e[i]]/n
    R = chol(sigma)  #cholesky decomp always exists for X'X/n where dj < n
    W.t = W[,pos_s[i]:pos_e[i]]%*%solve(R)
    return(list(R=R,W.t=W.t))
  }
  result <- lapply(seq(1,totalp),fun,W) 

  fun <- function(j,result)
  {
    result[[j]]$W.t
  }
  W.t <- matrix(unlist(lapply(seq(1,totalp),fun,result)),nrow=n)

  Y.hat <- X%*%solve(t(X)%*%X)%*%t(X)%*%Y
  r.alpha <- Y -Y.hat

  lambda.max = lammax(W.t,r.alpha,rep(1,totalp))$lambda.max
  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step )
  {
    lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  cv <- parCv(cl,lambda.seq,X,Y,W.t,group,fold.d,n.fold,gamma)
    
  list(lambda.seq=lambda.seq,PE=cv$PE,fold.d=fold.d)
}

#source("/space/jliu20/GE_interaction/cv.mcp.r")