sim <- function(num,p,q,mag1,mag2,mag,replic)
{
  tp <- numeric();tn <- numeric();fp <- numeric();fn <- numeric()
   
  for ( i in 1:replic)
  {
    data<-generate(mag,mag1,mag2,num,p,q,quan=0.9);
    X<-data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);

    index <- seq(1,(q+1)*p); #nonz <- index[data$b!=0];
    pvalue <- matrix(rep(0,(1+q)*p),nrow=(1+q))
    for ( j in 1:p)
    {
      fit <- lm(log(Y.s)~X.s+W.s[,index[(4*(j-1)+1):(4*j)]],weights=w.s)
      pvalue[,j] <- summary(fit)$coefficients[(2+q):(2+2*q),4]
    }
    pvalue <- as.vector(pvalue)
    s.gene <- index[pvalue < 5e-05]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/gene_s_match_surv_Marg_",num,"_",p,"_",mag1,"_",mag2,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);

    b <- data$b
    m.non <- match(s.gene,index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(s.gene)-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[pvalue >= 5e-05],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(s.gene)-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }

  list(tp=tp,fp=fp,tn=tn,fn=fn)
}

sim.band <- function(num,p,q,mag2,mag,replic)
{
  tp <- numeric();tn <- numeric();fp <- numeric();fn <- numeric()
 
  for ( i in 1:replic)
  {
    data <- generate(mag,mag2,num,p,q,quan=0.9);
    X<-data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);
  
    index <- seq(1,(q+1)*p); #nonz <- index[data$b!=0];
    pvalue <- matrix(rep(0,(1+q)*p),nrow=(1+q))
    for ( j in 1:p)
    {
      fit <- lm(log(Y.s)~X.s+W.s[,index[(4*(j-1)+1):(4*j)]],weights=w.s)
      pvalue[,j] <- summary(fit)$coefficients[(2+q):(2+2*q),4]
    }
    pvalue <- as.vector(pvalue)
    s.gene <- index[pvalue < 5e-05]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/gene_s_match_surv_Marg_",num,"_",p,"_",mag2,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);

    b <- data$b
    m.non <- match(s.gene,index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(s.gene)-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[pvalue >= 5e-05],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(s.gene)-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }

  list(tp=tp,fp=fp,tn=tn,fn=fn)
}

sim.band1 <- function(num,p,q,mag2,mag,replic)
{
  tp <- numeric();tn <- numeric();fp <- numeric();fn <- numeric()
    
  for ( i in 1:replic)
  {
    data <- generate(mag,mag2,num,p,q,quan=0.9);
    X<-data$Xs;Y<-data$Y;D<-data$D;W<-data$Ws;weight<-data$weight;T1 <- data$T1;
    w.s <- weight[weight!=0];X.s <- X[weight!=0,];W.s <-W[weight!=0,]; Y.s <-Y[weight!=0];n<-length(w.s);

    index <- seq(1,(q+1)*p); #nonz <- index[data$b!=0];
    pvalue <- matrix(rep(0,(1+q)*p),nrow=(1+q))
    for ( j in 1:p)
    {
      fit <- lm(log(Y.s)~X.s+W.s[,index[(4*(j-1)+1):(4*j)]],weights=w.s)
      pvalue[,j] <- summary(fit)$coefficients[(2+q):(2+2*q),4]
    }
    pvalue <- as.vector(pvalue)
    s.gene <- index[pvalue < 5e-05]
    outfile<- gsub("( )", "", paste("/space/jliu20/GE_interaction/Sim/gene_s_match_surv_b1_Marg_",num,"_",p,"_",mag2,".txt"))
    write.table(matrix(s.gene,nrow=1),outfile,sep=" ",quote=F, row.names=F,col.names=F,append=T);

    b <- data$b
    m.non <- match(s.gene,index[b!=0])
    tp <- c(tp,length(m.non[is.na(m.non)==F]) )
    fp <- c(fp,length(s.gene)-length(m.non[is.na(m.non)==F]) )
    m.zero <- match(index[pvalue >= 5e-05],index[b==0])
    tn <- c(tn, ncol(W) -length(b[b!=0])- (length(s.gene)-length(m.non[is.na(m.non)==F]))) 
    fn <- c(fn, length(b[b!=0])-length(m.non[is.na(m.non)==F]) )
  }

  list(tp=tp,fp=fp,tn=tn,fn=fn)
}
 
#source("/space/jliu20/GE_interaction/sim.trad.r")
#source("/space/jliu20/GE_interaction/generate.surv.r")
#s = sim(150,1000,3,0.2,0.05,1,5); #node01
#s = sim(150,1000,3,0.5,0.05,1,5); #node01
#s = sim(150,1000,3,0.8,0.05,1,5); #node01

#s = sim(250,1000,3,0.2,0.05,1,5); #node02
#s = sim(250,1000,3,0.5,0.05,1,5); #node02
#s = sim(250,1000,3,0.8,0.05,1,5); #node02


#source("/space/jliu20/GE_interaction/sim.trad.r")
#source("/space/jliu20/GE_interaction/generate.surv.band.r")
#s= sim.band(150,1000,3,0.05,1,5);  #node03
#s= sim.band(250,1000,3,0.05,1,5);  #node03

#source("/space/jliu20/GE_interaction/sim.trad.r")
#source("/space/jliu20/GE_interaction/generate.surv.band1.r")
#s= sim.band1(150,1000,3,0.05,1,5);   #node03
#s= sim.band1(250,1000,3,0.05,1,5);   #node03

