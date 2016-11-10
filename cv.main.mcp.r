.libPaths("~/Dropbox/workfolder/ABstudy/Codes-Liu")
library(snow)
home="~/Dropbox/workfolder/ABstudy/Codes-Liu/"

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

sse.evaluate <- function(i,X,Y,Z,fold.d,group,lambda,gamma)
{
  delete  <- unlist(fold.d[i])
  X.train <- X[-delete,]
  Y.train <- Y[-delete]
  Z.train <- Z[-delete,]
  
  X.valid <- X[delete,]
  Y.valid <- Y[delete]
  Z.valid <- Z[delete,]


  source(paste0(home, "CallC_mcp_main_iterate_wo.r"))
  fit <- mcp(X.train,Z.train,Y.train,lambda,gamma)
  
  sse <- sum((Y.valid - X.valid%*%fit$alpha - Z.valid%*%fit$beta)^2)
  #sse <- sum((y.valid -fit$beta0- x.valid%*%fit$beta)^2)
  list(sse=sse)
}

## looks like the main function 
parCv <- function(cl,lambda.seq,X,Y,Z,group,fold.d,n.fold,gamma)
{
  step <- length(lambda.seq)
  PE <- numeric(step)
  index <- seq(1:n.fold)

  for( i in 1:step)
  {
    lambda=lambda.seq[i]
    v <- lapply(index, sse.evaluate, X,Y,Z,fold.d,group,lambda,gamma)
    #v <- clusterApply(cl, index, sse.evaluate, X,Y,Z,fold.d,group,lambda,gamma)   #parallel
    sse.evaluate(i,X,Y,Z.t,fold.d,group,lambda,gamma)
    PE[i] <- sum(unlist(v))/n.fold
  }
  list(PE=PE)
}

#cv.optim <- function(X,Y,Z,group,n.fold,epsilon1,epsilon2,alpha,n.step,gamma,cl)
cv.optim <- function(X,Y,Z,n.fold,epsilon,n.step,gamma,cl){
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(Z)
  
  group <- rep(1,p)
  totalp <- sum(group)
  
  fold.d = splitList(sample(c(1:n),n),n.fold)    
  
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

  Y.hat <- X%*%solve(t(X)%*%X)%*%t(X)%*%Y
  r.alpha <- Y -Y.hat

  lambda.max = lammax(Z.t,r.alpha,rep(1,totalp))$lambda.max
  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step )
  {
    lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  cv <- parCv(cl,lambda.seq,X,Y,Z.t,group,fold.d,n.fold,gamma)
    
  list(lambda.seq=lambda.seq,PE=cv$PE,fold.d=fold.d)
}


