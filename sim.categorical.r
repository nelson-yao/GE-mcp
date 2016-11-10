
sim <- function(num,p,q,mag1,mag2,mag,gamma, epsilon, replic)
{
  tp <- numeric()
  tn <- numeric()
  fp <- numeric()
  fn <- numeric()
  lambda.opt.seq <- numeric()
  lambda.min.seq <- numeric()
  lambda.max.seq <- numeric()
  sse.seq <- numeric()
  
  for ( i in 1:replic)
  {
    data<-generate(mag,mag1,mag2,num,p,q,quan=0.9);
    X<-data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);

    weighted.standard.x <- function (x,w) {
       p <- ncol(x)
       n <- nrow(x)
       x.wmean <- matrix(rep(apply(x,2,weighted.mean,w),n),n,p,byrow=T)
       x.std <- (x-x.wmean)*w^0.5
    }
    weighted.standard.y <- function (y,w) {
       n <- length(y)
       y.wmean <- weighted.mean(y,w)
       y.std <- (y-y.wmean)*w^0.5
    }

    transform <- function(X,Y,W,weight)
    {
      X.sd=weighted.standard.x(X,weight); 
      W.sd=weighted.standard.x(W,weight);
      Y.sd=weighted.standard.y(Y,weight);  
      list(X.sd=X.sd,Y.sd=Y.sd,W.sd=W.sd)
    }

    data1=transform(X.s,log(Y.s),W.s,w.s);

    X=data1$X.sd; Y=data1$Y.sd; W=data1$W.sd;

    group <- rep(q+1,p)
    totalp <- sum(group)
    pos_e <- cumsum(group)
    pos_s <- pos_e - group + 1

    source("/space/jliu20/GE_interaction/cv1.r")

    n.fold <- 5; alpha <- 0.5; n.step <- 50; 
    cv <- cv.optim(X,Y,W,group,n.fold,epsilon,alpha,n.step,gamma,cl)
    lambda1.seq <- cv$lambda1.seq; lambda2.seq <- cv$lambda2.seq ;PE <- cv$PE
    lambda1.opt <- lambda1.seq[PE==min(PE)]; lambda2.opt <- lambda2.seq[PE==min(PE)]

    lambda.opt.seq <- c(lambda.opt.seq,lambda1.opt)
    lambda.min.seq <- c(lambda.min.seq,min(lambda1.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda1.seq))

    source("/space/jliu20/GE_interaction/CallC_iterate.r")
    fit <- GE_interact(Y,X,W,lambda1.opt,lambda2.opt,gamma,gamma)
    
    sse <- sum((Y-X%*%fit$alpha-fit$W.t%*%fit$beta)^2)
    sse.seq <- c(sse.seq, sse)  

    index=seq(1,ncol(W)); 
    s.gene=index[fit$beta!=0]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/gene_s_match_surv_sparse_gMCP_",num,"_",p,"_",mag1,"_",mag2,"_",gamma,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    
    b <- data$b
    m.non <- match(index[fit$beta!=0],index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[fit$beta==0],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }
  
  outfile1 <- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/sse_surv_sparse_gMCP_",num,"_",p,"_",mag1,"_",mag2,"_",gamma,".txt"))
  write.table(sse.seq,outfile1,sep=",",quote=F, row.names=F,col.names=F,append=T);

  list(tp=tp,fp=fp,tn=tn,fn=fn,sse.seq=sse.seq,lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}

sim.band <- function(num,p,q,mag2,mag,gamma, epsilon, replic)
{
  tp <- numeric()
  tn <- numeric()
  fp <- numeric()
  fn <- numeric()
  lambda.opt.seq <- numeric()
  lambda.min.seq <- numeric()
  lambda.max.seq <- numeric()
  sse.seq <- numeric()
  
  for ( i in 1:replic)
  {
    data <- generate(mag,mag2,num,p,q,quan=0.9);
    X <- data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);

    weighted.standard.x <- function (x,w) {
       p <- ncol(x)
       n <- nrow(x)
       x.wmean <- matrix(rep(apply(x,2,weighted.mean,w),n),n,p,byrow=T)
       x.std <- (x-x.wmean)*w^0.5
    }
    weighted.standard.y <- function (y,w) {
       n <- length(y)
       y.wmean <- weighted.mean(y,w)
       y.std <- (y-y.wmean)*w^0.5
    }

    transform <- function(X,Y,W,weight)
    {
      X.sd=weighted.standard.x(X,weight); 
      W.sd=weighted.standard.x(W,weight);
      Y.sd=weighted.standard.y(Y,weight);  
      list(X.sd=X.sd,Y.sd=Y.sd,W.sd=W.sd)
    }

    data1=transform(X.s,log(Y.s),W.s,w.s);

    X=data1$X.sd; Y=data1$Y.sd; W=data1$W.sd;

    group <- rep(q+1,p)
    totalp <- sum(group)
    pos_e <- cumsum(group)
    pos_s <- pos_e - group + 1

    #Sparse Group MCP
    source("/space/jliu20/GE_interaction/cv1.r")

    n.fold <- 5; alpha <- 0.5; n.step <- 50; 
    cv <- cv.optim(X,Y,W,group,n.fold,epsilon,alpha,n.step,gamma,cl)
    lambda1.seq <- cv$lambda1.seq; lambda2.seq <- cv$lambda2.seq ;PE <- cv$PE
    lambda1.opt <- lambda1.seq[PE==min(PE)]; lambda2.opt <- lambda2.seq[PE==min(PE)]

    lambda.opt.seq <- c(lambda.opt.seq,lambda1.opt)
    lambda.min.seq <- c(lambda.min.seq,min(lambda1.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda1.seq))

    source("/space/jliu20/GE_interaction/CallC_iterate.r")
    fit <- GE_interact(Y,X,W,lambda1.opt,lambda2.opt,gamma,gamma)
    
    sse <- sum((Y-X%*%fit$alpha-fit$W.t%*%fit$beta)^2)
    sse.seq <- c(sse.seq, sse)  

    index=seq(1,ncol(W)); 
    s.gene=index[fit$beta!=0]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/gene_s_match_surv_sparse_gMCP_",num,"_",p,"_",mag2,"_",gamma,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    
    b <- data$b
    m.non <- match(index[fit$beta!=0],index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[fit$beta==0],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }
  
  outfile1 <- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/sse_surv_sparse_gMCP_",num,"_",p,"_",mag2,"_",gamma,".txt"))
  write.table(sse.seq,outfile1,sep=",",quote=F, row.names=F,col.names=F,append=T);

  list(tp=tp,fp=fp,tn=tn,fn=fn,sse.seq=sse.seq,lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}

sim.band1 <- function(num,p,q,mag2,mag,gamma, epsilon, replic)
{
  tp <- numeric()
  tn <- numeric()
  fp <- numeric()
  fn <- numeric()
  lambda.opt.seq <- numeric()
  lambda.min.seq <- numeric()
  lambda.max.seq <- numeric()
  sse.seq <- numeric()
  
  for ( i in 1:replic)
  {
    data <- generate(mag,mag2,num,p,q,quan=0.9);
    X <- data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);

    weighted.standard.x <- function (x,w) {
       p <- ncol(x)
       n <- nrow(x)
       x.wmean <- matrix(rep(apply(x,2,weighted.mean,w),n),n,p,byrow=T)
       x.std <- (x-x.wmean)*w^0.5
    }
    weighted.standard.y <- function (y,w) {
       n <- length(y)
       y.wmean <- weighted.mean(y,w)
       y.std <- (y-y.wmean)*w^0.5
    }

    transform <- function(X,Y,W,weight)
    {
      X.sd=weighted.standard.x(X,weight); 
      W.sd=weighted.standard.x(W,weight);
      Y.sd=weighted.standard.y(Y,weight);  
      list(X.sd=X.sd,Y.sd=Y.sd,W.sd=W.sd)
    }

    data1=transform(X.s,log(Y.s),W.s,w.s);

    X=data1$X.sd; Y=data1$Y.sd; W=data1$W.sd;

    group <- rep(q+1,p)
    totalp <- sum(group)
    pos_e <- cumsum(group)
    pos_s <- pos_e - group + 1

    source("/space/jliu20/GE_interaction/cv1.r")

    n.fold <- 5; alpha <- 0.5; n.step <- 50; 
    cv <- cv.optim(X,Y,W,group,n.fold,epsilon,alpha,n.step,gamma,cl)
    lambda1.seq <- cv$lambda1.seq; lambda2.seq <- cv$lambda2.seq ;PE <- cv$PE
    lambda1.opt <- lambda1.seq[PE==min(PE)]; lambda2.opt <- lambda2.seq[PE==min(PE)]

    lambda.opt.seq <- c(lambda.opt.seq,lambda1.opt)
    lambda.min.seq <- c(lambda.min.seq,min(lambda1.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda1.seq))

    source("/space/jliu20/GE_interaction/CallC_iterate.r")
    fit <- GE_interact(Y,X,W,lambda1.opt,lambda2.opt,gamma,gamma)
    
    sse <- sum((Y-X%*%fit$alpha-fit$W.t%*%fit$beta)^2)
    sse.seq <- c(sse.seq, sse)  

    index=seq(1,ncol(W)); 
    s.gene=index[fit$beta!=0]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/gene_s_match_surv_b1_sparse_gMCP_",num,"_",p,"_",mag2,"_",gamma,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    
    b <- data$b
    m.non <- match(index[fit$beta!=0],index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[fit$beta==0],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(index[fit$beta!=0])-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }
  
  outfile1 <- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/Categorical/sse_surv_b1_sparse_gMCP_",num,"_",p,"_",mag2,"_",gamma,".txt"))
  write.table(sse.seq,outfile1,sep=",",quote=F, row.names=F,col.names=F,append=T);

  list(tp=tp,fp=fp,tn=tn,fn=fn,sse.seq=sse.seq,lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}

library(snow)
cl <- makeCluster(5)


#source("/space/jliu20/GE_interaction/sim.categorical.r")
#source("/space/jliu20/GE_interaction/generate.surv.categorical.r")
#s = sim(150,1000,3,0.2,0.05,1,6,0.02,5); #node06
#s = sim(150,1000,3,0.5,0.05,1,6,0.02,5); #node07
#s = sim(150,1000,3,0.8,0.05,1,6,0.01,5); #node08

#s = sim(250,1000,3,0.2,0.05,1,6,0.02,5); #node09
#s = sim(250,1000,3,0.5,0.05,1,6,0.02,5); #node11
#s = sim(250,1000,3,0.8,0.05,1,6,0.01,5); #node13

#source("/space/jliu20/GE_interaction/sim.categorical.r")
#source("/space/jliu20/GE_interaction/generate.surv.band.categorical.r")
#s= sim.band(150,1000,3,0.05,1,6,0.015,5); #node14
#s= sim.band(250,1000,3,0.05,1,6,0.01,5);  #node15

#source("/space/jliu20/GE_interaction/sim.categorical.r")
#source("/space/jliu20/GE_interaction/generate.surv.band1.categorical.r")
#s= sim.band1(150,1000,3,0.05,1,6,0.01,5);   #node2
#s= sim.band1(250,1000,3,0.05,1,6,0.01,5);   #node3

#node01 node02 node03 node04 node20
