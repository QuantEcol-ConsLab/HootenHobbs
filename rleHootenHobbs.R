# Import data
wt <- read.csv("C:/Users/Robbie/Desktop/HootenHobbsSupp/wt.csv")
setwd("C:/Users/Robbie/Desktop/HootenHobbsSupp")
library(jagsUI)

elev <- scale(wt$elev)
forest <- scale(wt$forest)

# Write model file
cat(file = "wtoccfull.txt",
    "model{
  # Priors
  b0 ~ dnorm(0, 0.444)
  b1 ~ dnorm(0, 0.444)
  b2 ~ dnorm(0, 0.444)
  p ~ dbeta(1,1)

  for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 + b1*elev[i] + b2*forest[i]
    }
    
    }")

# Create input data
M <- nrow(wt)
y <- rowSums(wt[,1:3],na.rm=T)
J <- apply(wt[,1:3], 1, function(x){length(which(!is.na(x)))})

win.data <- list(M = M, J = J, y = y, elev = elev, forest = forest)
#zinits <- ifelse(apply(y,1,sum)>0,1,0)
inits <- function(){list(p=runif(1))}
params <- c("b0", "b1", "b2", "p")
wt.out.full <- jagsUI::jags(win.data, inits, params, "wtoccfull.txt", n.chains=3, n.iter=5000,
               n.burnin=1000, parallel = TRUE, n.cores=3)

# Model with forest only
cat(file = "wtoccforest.txt",
    "model{
    # Priors
    b0 ~ dnorm(0, 0.444)
    b2 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 + b2*forest[i]
    }
    
    }")

params <- c("b0", "b2", "p")
wt.out.forest <- jagsUI::jags(win.data, inits, params, "wtoccforest.txt", n.chains=3, n.iter=5000,
                    n.burnin=1000, parallel = TRUE, n.cores=3)

# Model with elevation only
cat(file = "wtoccelev.txt",
    "model{
    # Priors
    b0 ~ dnorm(0, 0.444)
    b1 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0 + b1*elev[i]
    }
    
    }")

params <- c("b0", "b1", "p")
wt.out.elev <- jagsUI::jags(win.data, inits, params, "wtoccelev.txt", n.chains=3, n.iter=5000,
                      n.burnin=1000, parallel = TRUE, n.cores=3)

# Intercept only model

cat(file = "wtoccint.txt",
    "model{
    # Priors
    b0 ~ dnorm(0, 0.444)
    p ~ dbeta(1,1)
    
    for(i in 1:M){
    y[i]~dbin(p*z[i],J[i]) 
    
    z[i]=step(v[i])
    v[i]~dnorm(mu[i], 1)
    mu[i]=b0
    }
    
    }")

params <- c("b0", "p", "z")
wt.out.int <- jagsUI::jags(win.data, inits, params, "wtoccint.txt", n.chains=3, n.iter=5000,
                      n.burnin=1000, parallel = TRUE, n.cores=3)

## Bayes factors


