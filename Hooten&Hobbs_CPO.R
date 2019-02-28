## Hooten and Hobbs Model and Measures ##

setwd("~/Documents/Papers/Modeling/Hooten&Hobbs2015_Supplemental")

library(jagsUI);library(combinat)


## Prepare data ##

wt <- read.csv("wt.csv")
wt <- wt[1:200,]
## Simplify by summing across all surveys
y <- rowSums(wt[,1:3], na.rm = T)
elev <- as.vector(scale(wt$elev))
forest <- as.vector(scale(wt$forest))
M <- nrow(wt)
## Calculate number of survey occasions based on observations
J <- apply(wt[,1:3], 1, function(x){length(which(!is.na(x)))})

## Create empty matrix to populate with values as we find them
models <- c(1:4)
names <- c("Forest + Elev","Forest","Elev","Intercept")
cols <- c("-sum(logCPO)")
CPO_all <- matrix(NA, nrow = length(models), ncol = length(cols))
colnames(CPO_all) <- cols
rownames(CPO_all) <- names



### Set up occupancy model ###

### Full model ###

cat(file="HootinAndHobbinFULL.txt", 
"model {

  ## priors
  p ~ dbeta(1,1)
  b0 ~ dnorm(0,0.4444)
  b1 ~ dnorm(0,0.4444)
  b2 ~ dnorm(0,0.4444)
    
  for(i in 1:M){

    y[i] ~ dbin(p*z[i],J[i])

    # z[i] = step(v[i])
    # v[i] ~ dnorm(mu[i], 1)

    #Using phi() for the probit may be easier since we need psi[i] for CPO
    z[i] ~ dbin(psi[i], 1) 
    psi[i] <- phi(mu[i])         # probit link
    mu[i] = b0 + b1*elev[i] + b2*forest[i]

    ## Following H&H dmix function in occ.aux.mcmc.R script
    zero.idx[i] <- 1 - step(y[i]-0.1)   # 1 if no observations, 0 if observations
    ppd[i] <-  ((1-psi[i])+(psi[i]*(1-p)^J[i]))*zero.idx[i]  +  #if y==0, if y>0
               (1-zero.idx[i])*(psi[i]*
               (p^y[i] * (1-p)^(J[i]-y[i])* (exp(logfact(J[i]))/(exp(logfact(y[i]))*(exp(logfact(J[i]-y[i])))))))

    #monitor inverse ppd
    inv.ppd[i] <- 1/ppd[i]

  }
}",fill = TRUE)

## JAGS arguments and inits for N
data = list(M=M, J=J, elev=elev, forest=forest, y=y)
params = c('b0','b1','b2','p','inv.ppd') #'ppo.term','cpo.term')
inits = function() {list(p=runif(1), z=rep(1, M))}

nc = 3
ni = 2000
nb = 500
nthin = 1

out.full <- jagsUI::jags(data=data, inits=inits, model.file="HootinAndHobbinFULL.txt", parameters.to.save=params, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin)
#plot(out.full$samples[,"b0"])

#derive CPO[i] and -sum(log(CPO)) to match H&H Table 4
outmat <- as.matrix(out.full$samples)
inv.ppd <- outmat[,grep("inv.ppd", colnames(outmat))] #matrix of inv.ppd
CPO.vec <- (dim(outmat)[1])/apply(inv.ppd,2,sum)      #CPO[i]. H&H eq(21)
CPO_all[1,1] <- -sum(log(CPO.vec))                    #-log(CPO) in H&H Table 4

#Notice, this should be the similar to (identical to?)...
CPO.vec_v2 <-1/out.full$mean$inv.ppd       #posterior mean for each i
-sum(log(CPO.vec_v2))                    #should match above


save(out.full, file = "FullModel.Rdata")


### Forest only model ###

cat(file="HootinAndHobbinFOREST.txt", 
  "model {
    
  ## priors
  p ~ dbeta(1,1)
  b0 ~ dnorm(0,0.4444)
  b2 ~ dnorm(0,0.4444)
    
  for(i in 1:M){
    
    y[i] ~ dbin(p*z[i],J[i])
    
    # z[i] = step(v[i])
    # v[i] ~ dnorm(mu[i], 1)

    #Using phi() for the probit may be easier since we need psi[i] for CPO
    z[i] ~ dbin(psi[i], 1) 
    psi[i] <- phi(mu[i])         # probit link
    mu[i] = b0 + b2*forest[i]
  
    ## Following H&H dmix function in occ.aux.mcmc.R script
    zero.idx[i] <- 1 - step(y[i]-0.1)   # 1 if no observations, 0 if observations
    ppd[i] <-  ((1-psi[i])+(psi[i]*(1-p)^J[i]))*zero.idx[i]  +  #if y==0, if y>0
               (1-zero.idx[i])*(psi[i]* (p^y[i] * (1-p)^(J[i]-y[i])* (exp(logfact(J[i]))/(exp(logfact(y[i]))*(exp(logfact(J[i]-y[i])))))))
  
    #monitor inverse ppd
    inv.ppd[i] <- 1/ppd[i]
    
  }
}")

## JAGS arguments and inits for N
data = list(M=M, J=J, forest=forest, y=y)
params = c('b0','b2','p','inv.ppd')
inits = function() {list(p=runif(1), z=rep(1, M))}

nc = 3
ni = 2000
nb = 500
nthin = 1

out.forest <- jagsUI::jags(data=data, inits=inits, model.file="HootinAndHobbinFOREST.txt", parameters.to.save=params, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin)

#derive CPO[i] and -sum(log(CPO)) to match H&H Table 4
outmat <- as.matrix(out.forest$samples)
inv.ppd <- outmat[,grep("inv.ppd", colnames(outmat))] #matrix of inv.ppd
CPO.vec <- (dim(outmat)[1])/apply(inv.ppd,2,sum)      #CPO[i]. H&H eq(21)
CPO_all[2,1] <- -sum(log(CPO.vec))                    #-log(CPO) in H&H Table 4

#Notice, this should be the similar to (identical to?)...
CPO.vec_v2 <-1/out.forest$mean$inv.ppd       #posterior mean for each i
-sum(log(CPO.vec_v2))                    #should match above


save(out.forest, file = "ForestModel.Rdata")


## Elevation only model ##

cat(file="HootinAndHobbinELEV.txt", 
  "model {
    
  ## priors
  p ~ dbeta(1,1)
  b0 ~ dnorm(0,0.4444)
  b1 ~ dnorm(0,0.4444)
    
  for(i in 1:M){
    
    y[i] ~ dbin(p*z[i],J[i])
    
    # z[i] = step(v[i])
    # v[i] ~ dnorm(mu[i], 1)
  
    #Using phi() for the probit may be easier since we need psi[i] for CPO
    z[i] ~ dbin(psi[i], 1) 
    psi[i] <- phi(mu[i])         # probit link
    mu[i] = b0 + b1*elev[i]
  
    ## Following H&H dmix function in occ.aux.mcmc.R script
    zero.idx[i] <- 1 - step(y[i]-0.1)   # 1 if no observations, 0 if observations
    ppd[i] <-  ((1-psi[i])+(psi[i]*(1-p)^J[i]))*zero.idx[i]  +  #if y==0, if y>0
               (1-zero.idx[i])*(psi[i]* (p^y[i] * (1-p)^(J[i]-y[i])* (exp(logfact(J[i]))/(exp(logfact(y[i]))*(exp(logfact(J[i]-y[i])))))))
  
    #monitor inverse ppd
    inv.ppd[i] <- 1/ppd[i]
    
  }
}")

## JAGS arguments and inits for N
data = list(M=M, J=J, elev=elev, y=y)
params = c('b0','b1','p','inv.ppd')
inits = function() {list(p=runif(1), z=rep(1,M))}

nc = 3
ni = 2000
nb = 500
nthin = 1

out.elev <- jagsUI::jags(data=data, inits=inits, model.file="HootinAndHobbinELEV.txt", parameters.to.save=params, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin)

#derive CPO[i] and -sum(log(CPO)) to match H&H Table 4
outmat <- as.matrix(out.elev$samples)
inv.ppd <- outmat[,grep("inv.ppd", colnames(outmat))] #matrix of inv.ppd
CPO.vec <- (dim(outmat)[1])/apply(inv.ppd,2,sum)      #CPO[i]. H&H eq(21)
CPO_all[3,1] <- -sum(log(CPO.vec))                    #-log(CPO) in H&H Table 4

#Notice, this should be the similar to (identical to?)...
CPO.vec_v2 <-1/out.elev$mean$inv.ppd       #posterior mean for each i
-sum(log(CPO.vec_v2))                    #should match above

save(out.elev, file = "ElevModel.Rdata")


## Intercept only model ##

cat(file="HootinAndHobbinINT.txt", 
  "model {
    
  ## priors
  p ~ dbeta(1,1)
  b0 ~ dnorm(0,0.4444)
    
  for(i in 1:M){
    
    y[i] ~ dbin(p*z[i],J[i])
    
    # z[i] = step(v[i])
    # v[i] ~ dnorm(mu[i], 1)
  
    #Using phi() for the probit may be easier since we need psi[i] for CPO
    z[i] ~ dbin(psi[i], 1) 
    psi[i] <- phi(mu[i])         # probit link
    mu[i] = b0
  
    ## Following H&H dmix function in occ.aux.mcmc.R script
    zero.idx[i] <- 1 - step(y[i]-0.1)   # 1 if no observations, 0 if observations
    ppd[i] <-  ((1-psi[i])+(psi[i]*(1-p)^J[i]))*zero.idx[i]  +  #if y==0, if y>0
               (1-zero.idx[i])*(psi[i]* (p^y[i] * (1-p)^(J[i]-y[i])* (exp(logfact(J[i]))/(exp(logfact(y[i]))*(exp(logfact(J[i]-y[i])))))))
  
    #monitor inverse ppd
    inv.ppd[i] <- 1/ppd[i]
    
  }
}")

## JAGS arguments and inits for N
data = list(M=M, J=J, y=y)
params = c('b0','p', 'inv.ppd')
inits = function() {list(p=runif(1), z=rep(1,M))}

nc = 3
ni = 2000
nb = 500
nthin = 1

out.int <- jagsUI::jags(data=data, inits=inits, model.file="HootinAndHobbinINT.txt", parameters.to.save=params, n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nthin)

#derive CPO[i] and -sum(log(CPO)) to match H&H Table 4
outmat <- as.matrix(out.int$samples)
inv.ppd <- outmat[,grep("inv.ppd", colnames(outmat))] #matrix of inv.ppd
CPO.vec <- (dim(outmat)[1])/apply(inv.ppd,2,sum)      #CPO[i]. H&H eq(21)
CPO_all[4,1] <- -sum(log(CPO.vec))                    #-log(CPO) in H&H Table 4

#Notice, this should be the similar to (identical to?)...
CPO.vec_v2 <-1/out.int$mean$inv.ppd       #posterior mean for each i
-sum(log(CPO.vec_v2))                    #should match above

save(out.int, file = "IntModel.Rdata")


as.data.frame(CPO_all)


#### NOTES #####

### Conditional Predictive Ordinates ###
# probability (or density) of an observation without that data included in the model fitting
# Big CPO = very likely to include, small CPO = outlier
# Pro is that you can calc this as you fit your model by using a harmonic mean
# Con is that the harmonic mean might be unstable but there are program ways around this, is a within-sample measure that is likely to be overly "optimistic"
# Coded in to each model's run above
# All that exp(logfact())... is just dbinom(y[i],J[i],p)
# Could just do this outside of the model run as well and it would likely be faster and more elegant 

