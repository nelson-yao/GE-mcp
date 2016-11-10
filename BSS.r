source("/space/jliu20/GE_interaction/CallC_iterate.r")
   
sgm.bss <- function(X,Y,W,lambda1.opt,lambda2.opt,gamma,b,outfile){
   fit <- GE_interact(Y,X,W,lambda1.opt,lambda2.opt,gamma,gamma)
   index <- seq(1,ncol(W)); s.gene <- index[fit$beta!=0];

   n <- nrow(W); gene.s <- vector("list",b);
   
   for ( i in 1: b){
     sam <- sample(n,n*0.9); Y.sam <- Y[sam]; X.sam <- X[sam,]; W.sam <- W[sam,];
     fit.s <- GE_interact(Y.sam,X.sam,W.sam,lambda1.opt,lambda2.opt,gamma,gamma);
     s.gene.s <- index[fit.s$beta!=0];  gene.s[[i]] <- s.gene.s
     write.table(matrix(gene.s[[i]],nrow=1),outfile, sep=",",quote=F, row.names=F,col.names=F,append=T);
   }
   list(s.gene.s=s.gene.s,s.gene=s.gene)
}

source("/space/jliu20/GE_interaction/Group_MCP/CallC_iterate_gMCP.r")

gm.bss <- function(X,Y,W,lambda.opt,gamma,b,outfile){
   fit <- gmcp(Y,X,W,lambda.opt,gamma);
   index <- seq(1,ncol(W)); s.gene <- index[fit$beta!=0];

   n <- nrow(W); gene.s <- vector("list",b);
   
   for ( i in 1: b){
     sam <- sample(n,n*0.75); Y.sam <- Y[sam]; X.sam <- X[sam,]; W.sam <- W[sam,];
     fit.s <- gmcp(Y.sam,X.sam,W.sam,lambda.opt,gamma);
     s.gene.s <- index[fit.s$beta!=0];  gene.s[[i]] <- s.gene.s
     write.table(s.gene.s,outfile, sep=",",quote=F, row.names=F,col.names=F,append=T);
   }
   list(s.gene.s=s.gene.s)
}

source("/space/jliu20/GE_interaction/CallC_mcp_iterate.r")

m.bss <- function(X,Y,W,lambda.opt,gamma,b,outfile){
   fit <- mcp(X,W,Y,lambda.opt,gamma) ;
   index <- seq(1,ncol(W)); s.gene <- index[fit$beta!=0];

   n <- nrow(W); gene.s <- vector("list",b);
   
   for ( i in 1: b){
     sam <- sample(n,n*0.75); Y.sam <- Y[sam]; X.sam <- X[sam,]; W.sam <- W[sam,];
     fit.s <- mcp(Y.sam,X.sam,W.sam,lambda.opt,gamma);
     s.gene.s <- index[fit.s$beta!=0];  gene.s[[i]] <- s.gene.s
     write.table(s.gene.s,outfile, sep=",",quote=F, row.names=F,col.names=F,append=T);
   }
   list(s.gene.s=s.gene.s)
}
