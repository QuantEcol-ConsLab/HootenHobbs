
############ IMPORT LIBRARIES ###########
library(jagsUI)
library(rjags)
library(coda)

############ MODEL FILES ###########
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
    mu[i]=b0 + b1*elev[i]
    }}")
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
    mu[i]=b0+b2*forest[i]
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
    mu[i]=b0 + b1*elev[i] + b2*forest[i]
    }}")

############ SET UP DATA ###########
wt <- read.csv("/Users/abbybratt/Desktop/QCons_HootenHobbs/wt.csv")
wt <- wt[1:200, ] # subset to data used by HH
#scale data
elev <-as.numeric(scale(wt$elev))
forest <- as.numeric(scale(wt$forest))
# Create input data
M <- nrow(wt)
y <- rowSums(wt[,1:3],na.rm=T)
J <- apply(wt[,1:3], 1, function(x){length(which(!is.na(x)))})
win.data <- list(M = M, J = J, y = y, elev = elev, forest = forest)
inits <- function(){list(p=runif(1))}

############ SET UP MCMC ###########
n.burnin = 1000
n.iter = 5000

# likelihood function
dmix <- function(y,J,p,psi){
  out=rep(0,length(y))
  zero.idx=(y==0)
  out[zero.idx]=1-psi[zero.idx]+psi[zero.idx]*(1-p)^J[zero.idx]
  out[!zero.idx]=psi[!zero.idx]*dbinom(y[!zero.idx],J[!zero.idx],p)
  out
}

############ RUN MODELS ###########
# null
params <- c("b0", "p", "mu", "y", "z")
wt.out.int <- jagsUI::jags(win.data, inits, params, "wtocc_inter.txt", 
                           n.chains=3, n.iter=n.iter, n.burnin=n.burnin,
                           parallel = TRUE, n.cores=3)
wt.out.int$summary
summary(wt.out.int)
fit.int <- wt.out.int$samples
# elevation
params <- c("b0", "b1", "p", "mu", "y", "z")
wt.out.elev <- jagsUI::jags(win.data, inits, params, "wtocc_elev.txt", 
                            n.chains=3, n.iter=5000, n.burnin=1000,
                            parallel = TRUE, n.cores=3)
wt.out.elev$summary
summary(wt.out.elev)
fit.elev <- wt.out.elev$samples
# forest
params <- c("b0", "b2", "p", "mu", "y", "z")
wt.out.forest <- jagsUI::jags(win.data, inits, params, "wtocc_for.txt", 
                              n.chains=3, n.iter=5000, n.burnin=1000, 
                              parallel = TRUE, n.cores=3)
wt.out.forest$summary
summary(wt.out.forest)
fit.forest <- wt.out.forest$samples
# full
params <- c("b0", "b1", "b2", "p", "mu", "y", "z")
wt.out.full <- jagsUI::jags(win.data, inits, params, "wtocc_full.txt", 
                            n.chains=3, n.iter=5000, n.burnin=1000, 
                            parallel = TRUE, n.cores=3)
wt.out.full$summary
summary(wt.out.full)
fit.full <- wt.out.full$samples

############ COMPUTE DIC ###########
# null
# get ps and mus from each iteration and chain
p.int=c(fit.int[[1]][,"p"], fit.int[[2]][,"p"], fit.int[[3]][,"p"])
mu.int = rbind(fit.int[[1]][,3:202], fit.int[[2]][,3:202], fit.int[[3]][,3:202])
# mean across iterations
mean.mu.int = apply(mu.int, 2, mean)
ppd.int <- matrix(NA, nrow = (n.iter - n.burnin)*3, ncol = length(y))
for(i in 1:((n.iter - n.burnin)*3)) {
  ppd.int[i, ] <- dmix(y,J, mean(p.int), pnorm(mu.int[i, ])) # probit link
}
lppd.int= log(ppd.int)
# compute deviance evaluated at the posterior mean
D.hat.int=-2*sum(log(dmix(y,J,mean(p.int), pnorm(mean.mu.int))))
# compute posterior mean deviance
D.bar.int=mean(-2*apply(lppd.int,1,sum)) # eqn 38
pD.int=D.bar.int-D.hat.int # eqn 37
DIC.int=D.hat.int+2*pD.int # eqn 36
# print pD and DIC
c(pD.int, DIC.int)
# elev
p.elev=c(fit.elev[[1]][,"p"], fit.elev[[2]][,"p"], fit.elev[[3]][,"p"])
mu.elev = rbind(fit.elev[[1]][,4:203], fit.elev[[2]][,4:203], fit.elev[[3]][,4:203])
mean.mu.elev = apply(mu.elev, 2, mean)
ppd.elev <- matrix(NA, nrow = (n.iter - n.burnin)*3, ncol = length(y))
for(i in 1:((n.iter - n.burnin)*3)) {
  ppd.elev[i, ] <- dmix(y,J, mean(p.elev), pnorm(mu.elev[i, ]))
}
lppd.elev= log(ppd.elev)
D.hat.elev=-2*sum(log(dmix(y,J,mean(p.elev), pnorm(mean.mu.elev))))
D.bar.elev=mean(-2*apply(lppd.elev,1,sum))
pD.elev=D.bar.elev-D.hat.elev 
DIC.elev=D.hat.elev+2*pD.elev
c(pD.elev, DIC.elev)
# forest
p.forest=c(fit.forest[[1]][,"p"], fit.forest[[2]][,"p"], fit.forest[[3]][,"p"])
mu.forest = rbind(fit.forest[[1]][,4:203], fit.forest[[2]][,4:203], fit.forest[[3]][,4:203])
mean.mu.forest = apply(mu.forest, 2, mean)
ppd.forest <- matrix(NA, nrow = (n.iter - n.burnin)*3, ncol = length(y))
for(i in 1:((n.iter - n.burnin)*3)) {
  ppd.forest[i, ] <- dmix(y,J, mean(p.forest), pnorm(mu.forest[i, ]))
}
lppd.forest= log(ppd.forest)
D.hat.forest=-2*sum(log(dmix(y,J,mean(p.forest), pnorm(mean.mu.forest))))
D.bar.forest=mean(-2*apply(lppd.forest,1,sum))
pD.forest=D.bar.forest-D.hat.forest 
DIC.forest=D.hat.forest+2*pD.forest
c(pD.forest, DIC.forest)
# full
p.full=c(fit.full[[1]][,"p"], fit.full[[2]][,"p"], fit.full[[3]][,"p"])
mu.full = rbind(fit.full[[1]][,5:204], fit.full[[2]][,5:204], fit.full[[3]][,5:204])
mean.mu.full = apply(mu.full, 2, mean)
ppd.full <- matrix(NA, nrow = (n.iter - n.burnin)*3, ncol = length(y))
for(i in 1:((n.iter - n.burnin)*3)) {
  ppd.full[i, ] <- dmix(y,J, mean(p.full), pnorm(mu.full[i, ]))
}
lppd.full= log(ppd.full)
D.hat.full=-2*sum(log(dmix(y,J,mean(p.full), pnorm(mean.mu.full))))
D.bar.full=mean(-2*apply(lppd.full,1,sum))
pD.full=D.bar.full-D.hat.full 
DIC.full=D.hat.full+2*pD.full
c(pD.full, DIC.full)

############ COMPUTE WAIC ###########
# null
# -􏰇2 times the log point-wise predictive score
tmp.log=log(apply(ppd.int,2,mean)) 
tmp.sum=-2*sum(tmp.log)
pD.1.int=2*sum(tmp.log-apply(lppd.int,2,mean)) # eqn 43
pD.2.int=sum(apply(lppd.int,2,var)) # eq 44
WAIC.1.int=tmp.sum+2*pD.1.int
WAIC.2.int=tmp.sum+2*pD.2.int # this one preferred
c(pD.1.int, WAIC.1.int)
c(pD.2.int, WAIC.2.int)
# elev
tmp.log=log(apply(ppd.elev,2,mean))
tmp.sum=-2*sum(tmp.log)
pD.1.elev=2*sum(tmp.log-apply(lppd.elev,2,mean))
pD.2.elev=sum(apply(lppd.elev,2,var))
WAIC.1.elev=tmp.sum+2*pD.1.elev
WAIC.2.elev=tmp.sum+2*pD.2.elev
c(pD.1.elev, WAIC.1.elev)
c(pD.2.elev, WAIC.2.elev)
# forest
tmp.log=log(apply(ppd.forest,2,mean))
tmp.sum=-2*sum(tmp.log)
pD.1.forest=2*sum(tmp.log-apply(lppd.forest,2,mean))
pD.2.forest=sum(apply(lppd.forest,2,var))
WAIC.1.forest=tmp.sum+2*pD.1.forest
WAIC.2.forest=tmp.sum+2*pD.2.forest
c(pD.1.forest, WAIC.1.forest)
c(pD.2.forest, WAIC.2.forest)
# full
tmp.log=log(apply(ppd.full,2,mean))
tmp.sum=-2*sum(tmp.log)
pD.1.full=2*sum(tmp.log-apply(lppd.full,2,mean))
pD.2.full=sum(apply(lppd.full,2,var))
WAIC.1.full=tmp.sum+2*pD.1.full
WAIC.2.full=tmp.sum+2*pD.2.full
c(pD.1.full, WAIC.1.full)
c(pD.2.full, WAIC.2.full)

############ COMPUTE D ###########
n <- 2200# null
pz.chain1.int = fit.int[[1]][,403:602] * fit.int[[1]][,"p"]
pz.chain2.int =  fit.int[[2]][,403:602] * fit.int[[2]][,"p"]
pz.chain3.int =  fit.int[[3]][,403:602] * fit.int[[3]][,"p"]
pz.int = rbind(pz.chain1.int, pz.chain2.int, pz.chain3.int)
# some simulated ys
y.int = matrix(NA, nrow = (n.iter - n.burnin) * 3, ncol = n)
for (i in 1:((n.iter - n.burnin) * 3)) {
  y.int[i,] <- rbinom(n,J,pz.int[i, ])
}
# means
pz.mean.int = apply(pz.int,2,mean)
y.mean.int = apply(y.int,2,mean)
sum.1=sum((y-y.mean.int)^2) # sum of (data - simulated)^2
y.var.int=apply(y.int,2,var) # var of simulated
sum.2=sum(y.var.int) # sum of vars
D.int=sum.1+sum.2 # eqn 49
D.int
# elev
pz.chain1.elev = fit.elev[[1]][,404:603] * fit.elev[[1]][,"p"]
pz.chain2.elev =  fit.elev[[2]][,404:603] * fit.elev[[2]][,"p"]
pz.chain3.elev =  fit.elev[[3]][,404:603] * fit.elev[[3]][,"p"]
pz.elev = rbind(pz.chain1.elev, pz.chain2.elev, pz.chain3.elev)
y.elev = matrix(NA, nrow = (n.iter - n.burnin) * 3, ncol = n)
for (i in 1:((n.iter - n.burnin) * 3)) {
  y.elev[i,] <- rbinom(n,J,pz.elev[i, ])
}
pz.mean.elev = apply(pz.elev,2,mean)
y.mean.elev = apply(y.elev,2,mean)
sum.1=sum((y-y.mean.elev)^2)
y.var.elev=apply(y.elev,2,var)
sum.2=sum(y.var.elev)
D.elev=sum.1+sum.2
D.elev
# forest
pz.chain1.forest = fit.forest[[1]][,404:603] * fit.forest[[1]][,"p"]
pz.chain2.forest =  fit.forest[[2]][,404:603] * fit.forest[[2]][,"p"]
pz.chain3.forest =  fit.forest[[3]][,404:603] * fit.forest[[3]][,"p"]
pz.forest = rbind(pz.chain1.forest, pz.chain2.forest, pz.chain3.forest)
y.forest = matrix(NA, nrow = (n.iter - n.burnin) * 3, ncol = n)
for (i in 1:((n.iter - n.burnin) * 3)) {
  y.forest[i,] <- rbinom(n,J,pz.forest[i, ])
}
pz.mean.forest = apply(pz.forest,2,mean)
y.mean.forest = apply(y.forest,2,mean)
sum.1=sum((y-y.mean.forest)^2)
y.var.forest=apply(y.forest,2,var)
sum.2=sum(y.var.forest)
D.forest=sum.1+sum.2
D.forest
# full
pz.chain1.full = fit.full[[1]][,405:604] * fit.full[[1]][,"p"]
pz.chain2.full =  fit.full[[2]][,405:604] * fit.full[[2]][,"p"]
pz.chain3.full =  fit.full[[3]][,405:604] * fit.full[[3]][,"p"]
pz.full = rbind(pz.chain1.full, pz.chain2.full, pz.chain3.full)
y.full = matrix(NA, nrow = (n.iter - n.burnin) * 3, ncol = n)
for (i in 1:((n.iter - n.burnin) * 3)) {
  y.full[i,] <- rbinom(n,J,pz.full[i, ])
}
pz.mean.full = apply(pz.full,2,mean)
y.mean.full = apply(y.full,2,mean)
sum.1=sum((y-y.mean.full)^2)
y.var.full=apply(y.full,2,var)
sum.2=sum(y.var.full)
D.full=sum.1+sum.2
D.full

############ REPLICATE TABLE 5 ###########
c(DIC.int, DIC.elev, DIC.forest, DIC.full)
c(WAIC.2.int, WAIC.2.elev, WAIC.2.forest, WAIC.2.full)
c(D.int, D.elev, D.forest, D.full)
# all are a smidge off but rankings same