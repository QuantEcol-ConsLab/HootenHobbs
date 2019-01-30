####
####  Load Packages

library(rjags)

##Data are simulated and include a covariate x
library(rjags)
set.seed(11235)

poisim<-function(n=100, a0=0, a1=.8, a2=0, a3=-1){
x1<-runif(n, -3, 3)
x2<-runif(n, -3, 3)
x3<-runif(n, -3, 3)

lam=exp(a0 + a1*x1+ a2*x2+ a3*x3)
y<-rpois(n, lam)
dat=as.data.frame(cbind(y, x1, x2, x3))
return(dat)
}

pdat=poisim()

plot(pdat$x3, pdat$y)


sink("pois.txt")
cat("
model {

for (i in 1:n){
       y[i]~dpois(mu[i])
       log(mu[i])<- a0           
        }
  a0~dnorm(0, .01)
   }

", fill = TRUE)
sink()

data<-list(y=pdat$y, n=length(pdat$y))
parameters<-c("a0", "mu")
inits<-function(){list(a0=rnorm(1))}                                                         

#runs jags 
mod <- jags.model("pois.txt", data, inits, n.chains=3, n.adapt=100)
update(mod, 300)
fit <- coda.samples(mod,parameters,n.iter=5000)
DIC.jags=dic.samples(mod,n.iter=5000,type="pD")
DIC.jags

summary(fit)
plot(fit)


####
####  Calculate DIC in R 
####


lambda=rbind(fit[[1]][,-1], fit[[2]][,-1], fit[[3]][,-1])
mean.lambda = apply(lambda,2,mean)
Dhat = -2*(sum(dpois(pdat$y,mean.lambda,log=TRUE))) 
Dbar = mean(-2*apply(dpois(pdat$y,t(lambda),log=TRUE),2,sum))
pD.DIC = Dbar - Dhat
DIC = Dhat + 2*pD.DIC
c(pD.DIC,DIC)

####
####  Calculate WAIC from JAGS Output 
####

lppd=-2*sum(log(apply(dpois(pdat$y,t(lambda)),1,mean)))
pD.WAIC=sum(apply(dpois(pdat$y,t(lambda),log=TRUE),1,var))
WAIC=lppd+2*pD.WAIC
c(pD.WAIC,WAIC)

####
####  Calculate D_inf from JAGS Output 
####

y.pred=matrix(0,length(pdat$y),15000)
for(k in 1:15000){
  y.pred[,k]=rpois(length(pdat$y),lambda[k,])
}  

D.inf=sum((pdat$y-apply(y.pred,1,mean))^2)+sum(apply(y.pred,1,var))
D.inf




####
####  Load Packages

library(rjags)

##Data are simulated and include a covariate x
library(rjags)
set.seed(11235)

poisim<-function(n=100, a0=0, a1=.8, a2=0, a3=-1){
x1<-runif(n, -3, 3)
x2<-runif(n, -3, 3)
x3<-runif(n, -3, 3)

lam=exp(a0 + a1*x1+ a2*x2+ a3*x3)
y<-rpois(n, lam)
dat=as.data.frame(cbind(y, x1, x2, x3))
return(dat)
}

pdat=poisim()

plot(pdat$x3, pdat$y)


sink("pois.txt")
cat("
model {

for (i in 1:n){
       y[i]~dpois(mu[i])
       log(mu[i])<- a0 +a1*x1[i]          
        }
  a0~dnorm(0, .01)
  a1~dnorm(0, .01)
   }

", fill = TRUE)
sink()

data<-list(y=pdat$y, x1=pdat$x1, n=length(pdat$y))
parameters<-c("a0", "a1", "mu")
inits<-function(){list(a0=rnorm(1), a1=rnorm(1))}                                                         

#runs jags 
mod <- jags.model("pois.txt", data, inits, n.chains=3, n.adapt=100)
update(mod, 300)
fit <- coda.samples(mod,parameters,n.iter=5000)

summary(fit)
plot(fit)




####
####  Calculate DIC from JAGS Output 
####

library(coda)


lambda=rbind(fit[[1]][,-(1:2)], fit[[2]][,-(1:2)], fit[[3]][,-(1:2)])
mean.lambda = apply(lambda,2,mean)
Dhat = -2*(sum(dpois(pdat$y,mean.lambda,log=TRUE))) 
Dbar = mean(-2*apply(dpois(pdat$y,t(lambda),log=TRUE),2,sum))
pD.DIC = Dbar - Dhat
DIC = Dhat + 2*pD.DIC
c(pD.DIC,DIC)

####
####  Calculate WAIC from JAGS Output 
####

lppd=-2*sum(log(apply(dpois(pdat$y,t(lambda)),1,mean)))
pD.WAIC=sum(apply(dpois(pdat$y,t(lambda),log=TRUE),1,var))
WAIC=lppd+2*pD.WAIC
c(pD.WAIC,WAIC)

####
####  Calculate D_inf from JAGS Output 
####

y.pred=matrix(0,length(pdat$y),15000)
for(k in 1:15000){
  y.pred[,k]=rpois(length(pdat$y),lambda[k,])
}  

D.inf=sum((pdat$y-apply(y.pred,1,mean))^2)+sum(apply(y.pred,1,var))
D.inf

####
#### Calculate DIC using Built-in JAGS Function 
####
#### This needs multiple chains to work.
####

DIC.jags=dic.samples(mod,n.iter=5000,type="pD")
DIC.jags


