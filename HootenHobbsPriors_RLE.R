#### Shrinkage from different priors ####
wt <- read.csv("C:/Users/Robbie/Desktop/HootenHobbsSupp/wt.csv")
setwd("C:/Users/Robbie/Desktop/wolverines/wolverinepractice/old_and_backup/")
library(jagsUI)
#scale data
elevation <-as.numeric(scale(wt$elev))
forestdata <- as.numeric(scale(wt$forest))
####models with the original prior #######
#Intercept only model
cat(file = "wtocc_inter.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
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
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
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
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
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
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i] + b2*forestdata[i]
    }}")

#### Models with vague normal priors on b ####
#Intercept only model
cat(file = "wtocc_inter_vague.txt",
    "model{
    b0 ~ dnorm(0, 0.01)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 
    }}")
#elevation only model
cat(file = "wtocc_elev_vague.txt",
    "model{
    b0 ~ dnorm(0, 0.01)
    b1 ~ dnorm(0, 0.01)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i]
    } }")
#forest only model
cat(file = "wtocc_for_vague.txt",
    "model{
    b0 ~ dnorm(0, 0.01)
    b2 ~ dnorm(0, 0.01)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0+b2*forestdata[i]
    }}")
#full model
cat(file = "wtocc_full_vague.txt",
    "model{
    b0 ~ dnorm(0, 0.01)
    b1 ~ dnorm(0, 0.01)
    b2 ~ dnorm(0, 0.01)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i] + b2*forestdata[i]
    }}")

#### Models with uniform priors on b ####
#Intercept only model
cat(file = "wtocc_inter_unif.txt",
    "model{
    b0 ~ dunif(-10, 10)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 
    }}")
#elevation only model
cat(file = "wtocc_elev_unif.txt",
    "model{
    b0 ~ dunif(-10, 10)
    b1 ~ dunif(-10, 10)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i]
    } }")
#forest only model
cat(file = "wtocc_for_unif.txt",
    "model{
    b0 ~ dunif(-10, 10)
    b2 ~ dunif(-10, 10)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0+b2*forestdata[i]
    }}")
#full model
cat(file = "wtocc_full_unif.txt",
    "model{
    b0 ~ dunif(-10, 10)
    b1 ~ dunif(-10, 10)
    b2 ~ dunif(-10, 10)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i] + b2*forestdata[i]
    }}")

#### Models with Laplace priors on b - like lasso ####
#Intercept only model
cat(file = "wtocc_inter_laplace.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 
    }}")
#elevation only model
cat(file = "wtocc_elev_laplace.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    tau ~ dunif(0.0001, 10)
    b1 ~ ddexp(0, tau)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i]
    } }")
#forest only model
cat(file = "wtocc_for_laplace.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    tau ~ dunif(0.0001, 10)
    b2 ~ ddexp(0, tau)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0+b2*forestdata[i]
    }}")
#full model
cat(file = "wtocc_full_laplace.txt",
    "model{
    b0 ~ dnorm(0, 0.444)
    tau ~ dunif(0.001, 10)
    b1 ~ ddexp(0, tau)
    b2 ~ ddexp(0, tau)
    p ~ dbeta(1,1)
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i] ~ dbin(psi[i], 1)
    psi[i] = phi(mu[i])
    mu[i]=b0 + b1*elevation[i] + b2*forestdata[i]
    }}")

#### Fit and compare models with different priors ####
## "normal" normal
M <- nrow(wt)
y <- rowSums(wt[,1:3],na.rm=T)
J <- apply(wt[,1:3], 1, function(x){length(which(!is.na(x)))})

win.data <- list(M = M, J = J, y = y, elevation = elevation, forestdata = forestdata)
#zinits <- ifelse(apply(y,1,sum)>0,1,0)
inits <- function(){list(p=runif(1), z=rep(1,M))}

params <- c("b0", "p")
wt.out.inter <- jagsUI::jags(win.data, inits, params, "wtocc_inter.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "p")
wt.out.elev <- jagsUI::jags(win.data, inits, params, "wtocc_elev.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b2", "p")
wt.out.for <- jagsUI::jags(win.data, inits, params, "wtocc_for.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "b2", "p")
wt.out.full <- jagsUI::jags(win.data, inits, params, "wtocc_full.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

## vague normal
params <- c("b0", "p")
wt.out.inter.vague <- jagsUI::jags(win.data, inits, params, "wtocc_inter_vague.txt", n.chains=3, n.iter=5000,
                             n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "p")
wt.out.elev.vague <- jagsUI::jags(win.data, inits, params, "wtocc_elev_vague.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b2", "p")
wt.out.for.vague <- jagsUI::jags(win.data, inits, params, "wtocc_for_vague.txt", n.chains=3, n.iter=5000,
                           n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "b2", "p")
wt.out.full.vague <- jagsUI::jags(win.data, inits, params, "wtocc_full_vague.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

## uniform
params <- c("b0", "p")
wt.out.inter.unif <- jagsUI::jags(win.data, inits, params, "wtocc_inter_unif.txt", n.chains=3, n.iter=5000,
                             n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "p")
wt.out.elev.unif <- jagsUI::jags(win.data, inits, params, "wtocc_elev_unif.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b2", "p")
wt.out.for.unif <- jagsUI::jags(win.data, inits, params, "wtocc_for_unif.txt", n.chains=3, n.iter=5000,
                           n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "b2", "p")
wt.out.full.unif <- jagsUI::jags(win.data, inits, params, "wtocc_full_unif.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

## Laplace
params <- c("b0", "p")
wt.out.inter.lap <- jagsUI::jags(win.data, inits, params, "wtocc_inter_laplace.txt", n.chains=3, n.iter=5000,
                             n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "p")
wt.out.elev.lap <- jagsUI::jags(win.data, inits, params, "wtocc_elev_laplace.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b2", "p")
wt.out.for.lap <- jagsUI::jags(win.data, inits, params, "wtocc_for_laplace.txt", n.chains=3, n.iter=5000,
                           n.burnin=1000, parallel = TRUE, n.cores=3)

params <- c("b0", "b1", "b2", "p")
wt.out.full.lap <- jagsUI::jags(win.data, inits, params, "wtocc_full_laplace.txt", n.chains=3, n.iter=5000,
                            n.burnin=1000, parallel = TRUE, n.cores=3)

## Comparison
cbind(wt.out.inter$mean, wt.out.inter.vague$mean, wt.out.inter.unif$mean, wt.out.inter.lap$mean)
cbind(wt.out.elev$mean, wt.out.elev.vague$mean, wt.out.elev.unif$mean, wt.out.elev.lap$mean)
cbind(wt.out.for$mean, wt.out.for.vague$mean, wt.out.for.unif$mean, wt.out.for.lap$mean)
cbind(wt.out.full$mean, wt.out.full.vague$mean, wt.out.full.unif$mean, wt.out.full.lap$mean)
