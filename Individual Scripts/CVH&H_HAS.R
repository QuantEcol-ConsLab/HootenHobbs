wt <- read.csv("/Users/sipeh/Desktop/Lab Meetings/wt.csv")
setwd("/Users/sipeh/Desktop/Lab Meetings")
library(jagsUI)
#scale data
elevation <-as.numeric(scale(wt$elev))
forestdata <- as.numeric(scale(wt$forest))
####models for functions#######
#Intercept only model
cat(file = "wtocc_inter.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 
    }}")
#elevation only model
cat(file = "wtocc_elev.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    b1 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 + b1*elevation[i]
    } }")
#forest only model
cat(file = "wtocc_for.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    b2 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0+b2*forestdata[i]
    }}")
#full model
cat(file = "wtocc_full.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    b1 ~ dnorm(0, 0.444)
    b2 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 + b1*elevation[i] + b2*forestdata[i]
    }}")


#function to run JAGs model and compute the CV score
#description of each input is described below
Jags_and_CV<-function(n,K, model, n.chains, n.iter, parallel, n.cores, n.burnin, params, inter, elev, forest, full){
  y <- rowSums(wt[,1:3],na.rm=T)
  J <- apply(wt[,1:3], 1, function(x){length(which(!is.na(x)))})
  fold_mat<-matrix(1:n, nrow=(n/K), ncol=K)
  data<-vector(mode="list", length=K)
  outs<-vector(mode="list", length=K)
  meanppd<-matrix(nrow=length(fold_mat[,1]), ncol=K)
  #choose which data to use in the model
  for(k in 1:K){
    if(inter==TRUE){
    data[[k]] <- list(M = (n-(n/K)), J = J[-fold_mat[,k]], y = y[-fold_mat[,k]])
    }else{}
    if(elev==TRUE){
      data[[k]] <- list(M = (n-(n/K)), J = J[-fold_mat[,k]], y = y[-fold_mat[,k]], 
                        elevation=elevation[-fold_mat[,k]])
    }else{}
    if(forest==TRUE){
      data[[k]] <- list(M = (n-(n/K)), J = J[-fold_mat[,k]], y = y[-fold_mat[,k]], 
                        forestdata=forestdata[-fold_mat[,k]])
    }else{}
    if(full==TRUE){
      data[[k]] <- list(M = (n-(n/K)), J = J[-fold_mat[,k]], y = y[-fold_mat[,k]], 
                        elevation=elevation[-fold_mat[,k]], forestdata=forestdata[-fold_mat[,k]])
    }else{}
    inits<-function(){list(p=runif(1))}
    #run JAGs model 
    outs[[k]] <- jags(data[[k]], inits, params, model, n.chains=n.chains, n.iter=n.iter,
                      n.burnin=n.burnin, parallel = parallel, n.cores=n.cores)
  }
  #compute the score for each model
  for(k in 1:K){
     if(inter==TRUE){
     meanppd[,k]<-dmix(y=y[fold_mat[,k]], J=J[fold_mat[,k]], p=outs[[k]]$mean$p,
                       psi=pnorm(rep(outs[[k]]$mean$b0, (n/K))))
      }else{}
    if(elev==TRUE){
      meanppd[,k]<-dmix(y=y[fold_mat[,k]], J=J[fold_mat[,k]], p=outs[[k]]$mean$p,
                        psi=pnorm(rep(outs[[k]]$mean$b0, (n/K))+elevation[fold_mat[,k]]*outs[[k]]$mean$b1))
    }else{}
    if(forest==TRUE){
      meanppd[,k]<-dmix(y=y[fold_mat[,k]], J=J[fold_mat[,k]], p=outs[[k]]$mean$p,
                        psi=pnorm(rep(outs[[k]]$mean$b0, (n/K))+ forestdata[fold_mat[,k]]*outs[[k]]$mean$b2))
    }else{}
    if(full==TRUE){
      meanppd[,k]<-dmix(y=y[fold_mat[,k]], J=J[fold_mat[,k]], p=outs[[k]]$mean$p,
                        psi=pnorm(rep(outs[[k]]$mean$b0, (n/K))+elevation[fold_mat[,k]]*outs[[k]]$mean$b1+
                        forestdata[fold_mat[,k]]*outs[[k]]$mean$b2))
    }else{}
  }
  score<- -2*sum(log(meanppd))
  
  return(list(score=score, meanppd=meanppd))#, outs=outs))
}
#function from H&H updated R code to compute the score function
dmix <- function(y,J,p,psi){
  out=rep(0,length(y))
  zero.idx=(y==0)
  out[zero.idx]=1-psi[zero.idx]+psi[zero.idx]*(1-p)^J[zero.idx]
  out[!zero.idx]=psi[!zero.idx]*dbinom(y[!zero.idx],J[!zero.idx],p)
  out
}

#definition for input to function:
#n is the total data, H&H used 200 so I did as well
#K is the number of folds, H&H used 10
#the way I wrote the function above, n/K must be a whole number
# choose the model you want to run by True or False
# n.chains, n.iter, n.burnin, n.cores, parallel (T/F) are all for running the model in JAGs

#run each model and compare scores
inter_model<-Jags_and_CV(n=200, K=10, model="wtocc_inter.txt", n.chains=3, n.iter=10000,
                          parallel=TRUE, n.cores=3, n.burnin=5000, params<-c("p","b0"),
                          inter=TRUE, elev=FALSE, forest=FALSE, full=FALSE) 
elev_model<-Jags_and_CV(n=200, K=10, model="wtocc_elev.txt", n.chains=3, n.iter=10000,
                          parallel=TRUE, n.cores=3, n.burnin=5000, params<-c("p","b0","b1"),
                          inter=FALSE, elev=TRUE, forest=FALSE, full=FALSE) 
forest_model<-Jags_and_CV(n=200, K=10, model="wtocc_for.txt", n.chains=3, n.iter=10000,
                        parallel=TRUE, n.cores=3, n.burnin=5000, params<-c("p","b0","b2"),
                        inter=FALSE, elev=FALSE, forest=TRUE, full=FALSE) 
full_model<-Jags_and_CV(n=200, K=10, model="wtocc_full.txt", n.chains=3, n.iter=10000,
                          parallel=TRUE, n.cores=3, n.burnin=5000, params<-c("p","b0","b1","b2"),
                          inter=FALSE, elev=FALSE, forest=FALSE, full=TRUE) 

#output results into a dataframe
results<-matrix(nrow=4, ncol=2)
colnames(results)<-c( "Covariates", "C-V Score")
rownames(results)<-c("M1", "M2", "M3", "M4")
results[,1]<-c("Null", "Elev", "For", "Full")
results[,2]<-c(signif(inter_model$score,4), signif(elev_model$score,4),
               signif(forest_model$score,4), signif(full_model$score,4))
as.data.frame(results)

