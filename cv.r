#.libPaths("~/rpacks")
library(snow)

lambda2.max <- function(x,y,group)
{
  n <- nrow(x)
  remove.id <- seq(1,sum(group),group[1])
  tmp <- abs(colSums(as.vector(y)*x))[-remove.id]
  lambda.max <- max(tmp)/n
  list(lambda2.max=lambda.max)
}

lambda1.max <- function(x,y,group,pos_s,pos_e,lambda2)
{
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
  p <- length(group)
  fun <- function(j,pos_s,pos_e,x,y,group,lambda2)
  {
    t=x[,pos_s[j]:pos_e[j]]
    n=nrow(t)
    z=colSums(t*as.vector(y)/n)
    tmp=(abs(z)-lambda2)
    tmp[1]=z[1]
    tmp[tmp<=0]=0
    s_value=sign(z)*tmp
    inner = sum(s_value^2)^0.5/group[j]^0.5
    #inner=sum( abs(colSums(t*as.vector(y)/n)) ^2)^0.5/group[j]^0.5
    list(inner=inner)
  }
  lambda.max=max(unlist(lapply(seq(1,p),fun,pos_s,pos_e,x,y,group,lambda2)))
  list(lambda1.max=lambda.max)
}

sse.evaluate <- function(i,X,Y,W,fold.d,group,lambda1,lambda2,gamma)
{
  delete  <- unlist(fold.d[i])
  X.train <- X[-delete,]
  Y.train <- Y[-delete]
  W.train <- W[-delete,]
  
  X.valid <- X[delete,]
  Y.valid <- Y[delete]
  W.valid <- W[delete,]

  source("/space/jliu20/GE_interaction/CallC_iterate.r")
  fit <- GE_interact(Y.train,X.train,W.train,lambda1,lambda2,gamma,gamma)

  p <- length(group)
  n <- nrow(W.valid)
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
  
  fun <- function(j,result,W,pos_s, pos_e)
  {
    R <- result[[j]]$R
    W[,pos_s[i]:pos_e[i]]%*%solve(R)
  }
  W.v <- matrix(unlist(lapply(seq(1,p),fun,fit$result, W.valid, pos_s,pos_e)),nrow=n)

  sse <- sum((Y.valid - X.valid%*%fit$alpha - W.v%*%fit$beta)^2)
  #sse <- sum((y.valid -fit$beta0- x.valid%*%fit$beta)^2)
  list(sse=sse)
}

parCv <- function(cl,lambda1.seq,lambda2.seq,X,Y,W,group,fold.d,n.fold,gamma)
{
  step <- length(lambda1.seq)
  PE <- numeric(step)
  index <- seq(1:n.fold)

  for( i in 1:step)
  {
    lambda1=lambda1.seq[i]
    lambda2=lambda2.seq[i]
    #v <- lapply(index,fun,X,Y,W,fold.d,n.fold,group,lambda1,lambda2,gamma)
    v <- clusterApply(cl, index, sse.evaluate, X,Y,W,fold.d,group,lambda1,lambda2,gamma)
    PE[i] <- sum(unlist(v))/n.fold
  }
  list(PE=PE)
}

#cv.optim <- function(X,Y,W,group,n.fold,epsilon1,epsilon2,alpha,n.step,gamma,cl)
cv.optim <- function(X,Y,W,group,n.fold,epsilon,alpha,n.step,gamma,cl)
{
  p <- length(group)
  n <- dim(W)[1]

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
  result <- lapply(seq(1,p),fun,W) 

  fun <- function(j,result)
  {
    result[[j]]$W.t
  }
  W.t <- matrix(unlist(lapply(seq(1,p),fun,result)),nrow=n)

  Y.hat <- X%*%solve(t(X)%*%X)%*%t(X)%*%Y
  r.alpha <- Y -Y.hat

  #l2.max <- lambda2.max(W.t,r.alpha,group)$lambda2.max
  #l2.min <- l2.max*epsilon2
  #ss2 <- (log(l2.max)-log(l2.min))/(n.step-1)
  
  #lambda1.seq <- numeric(n.step*n.step)
  #lambda2.seq <- numeric(n.step*n.step)
  
  #for ( i2 in 1:n.step )
  #{ 
  #  lambda2 <- exp(log(l2.max)-ss2*(i2-1)) 
  #  l1.max <- lambda1.max(W.t,r.alpha,group,pos_s,pos_e,lambda2)$lambda1.max
  #  l1.min <- l1.max*epsilon1
  #  ss1 <- (log(l1.max)-log(l1.min))/(n.step-1)
  #  for ( i1 in 1:n.step )
  #  {
  #     lambda1.seq[i1+(i2-1)*n.step] <- exp(log(l1.max)-ss1*(i1-1))
  #     lambda2.seq[i1+(i2-1)*n.step] <- lambda2
  #  }
  #}
 
  l2.max <- lambda2.max(W.t,r.alpha,group)$lambda2.max
  l1.max <- lambda1.max(W.t,r.alpha,group,pos_s,pos_e,0)$lambda1.max *1.5
  if ( alpha != 1)
  {
   lambda1.max <- l1.max
   lambda1.min <- lambda1.max*epsilon
   ss <- (log(lambda1.max)-log(lambda1.min))/(n.step-1)
   lambda2.seq <- numeric(n.step)
   lambda1.seq <- numeric(n.step)
   for ( i in 1:n.step )
   {
    lambda1.seq[i] <- exp(log(lambda1.max)-ss*(i-1))
    lambda2.seq[i] <- (1-alpha)/alpha*lambda1.seq[i]
   }
  }
  else if (alpha == 1)
  {
   lambda1.max <- l1.max
   lambda1.min <- lambda1.max*epsilon
   ss <- (log(lambda1.max)-log(lambda1.min))/(n.step-1)
   lambda2.seq <- numeric(n.step)
   lambda1.seq <- numeric(n.step)
   for ( i in 1:n.step )
   {
    lambda2.seq[i] <- 0
    lambda1.seq[i] <- exp(log(lambda1.max)-ss*(i-1))
   }
  }

  cv <- parCv(cl,lambda1.seq,lambda2.seq,X,Y,W,group,fold.d,n.fold,gamma)
    
  list(lambda1.seq=lambda1.seq,lambda2.seq=lambda2.seq,PE=cv$PE,fold.d=fold.d)
}

#source("/space/jliu20/GE_interaction/cv.r")