# Penalized-Lq for JLSSM-SN (one dimension)
rm(list=ls())
library(gdata) # for lower triangle
library(MASS)
library(mvtnorm)
library(DiceKriging) #SCAD penalty
library(sn)
library(lmtest)
options(scipen = 999) #disabling scientific notation 

#Generate data
n = 100  #number of observations

a = 1 #iteration number
b=1
c=1
gvec = NULL
g =0.5 # NULL  #Bridge search interval (0.01, 1.8) by =0.01
eps = 0.05 # epsilon
beta_diff = 10000 #to get in the loop at the beginning
gamma_diff = 10000 
alpha_diff = 10000 
maxiter = 100  #If the convergence is not achieved, it tells after how many iterations it will stop.

### Real parameters
beta0 = matrix(c(1, 1, 0, 0, 1, 0, 0), ncol=1)
gamma0 = matrix(c(0.7, 0.7, 0, 0, 0.7), ncol=1)
alpha0 = matrix(c(0.5, 0.5, 0, 0, 0.5), ncol=1)

p = length(beta0) #dimension of location model parameter vector
d = length(gamma0) #dimension of scale model parameter vector
q= length(alpha0) #dimension of skewness model parameter vector
d2=d+1
d3=d+q

### Initial values
beta=beta0
gamma=gamma0
alpha=alpha0

### Grid search for penalty terms 
penalty = 2 # 1 for LASSO ve BRIDGE, 2 for SCAD
#tau = seq(0.0001, 2, by=0.01) # 
tau = seq(0, 0.9, by=0.1) # search (0.1, 3) by=0.01 
#Lq=seq(0.01,1,by=0.01) interval search Lq
Lq=0.9 
BIC = matrix(NA, length(tau), 1)
result_BRDG <- NULL

### Weight matrix for penalty
W_beta <- matrix(0, p, p)
W_gamma <- matrix(0, d, d)
W_alpha <- matrix(0, q, q)

### Covariates 
X = matrix(1, n, p)
Z = matrix(1, n, d)
V = matrix(1, n, q)

for(j in 1:p){X[,j] <- matrix(rnorm(n, 0,1),ncol=1)}
for(j in 1:d){Z[,j] <- matrix(rnorm(n, 0,1),ncol=1)}
for(j in 1:q){V[,j] <- matrix(rnorm(n,0,1),ncol=1)}

### Distribution parameters and kappa and delta for EM
mu_i=X%*%beta0
sigma_i=exp((Z%*%gamma0)/2)
lambda_i=V%*%alpha0

delta_i=lambda_i/sqrt(1+(lambda_i)^2)
kappa_i2=((sigma_i)^2)*(1-(delta_i)^2)


# Replication Starts (Generate data - response)
y<-sapply(1:n, function(i) rsn(1,mu_i[i],sigma_i[i],lambda_i[i])) # y_i <- y[,i]
#y <- sapply(1:n, function(i) mu_i[i]+sigma_i[i]*delta_i[i]*abs(rnorm(1))+sigma_i[i]*sqrt(1-delta_i[i]^2)*rnorm(1)) # y_i <- y[,i]
r <- y-mu_i  #residuals

# or 
#lm.obj <- lm(y ~ X - 1)
#beta <- coef(lm.obj)
#resid(lm.obj) -> res
#gamma <- coef(lm(log(res^2) ~ Z - 1))
#alpha <- rep(0, q)

####################################
jlssm <- function(beta, gamma, alpha, X, V, Z, y, Lq, penalty, g,  tt){
  #cat("firstbeta", beta, "\n")
  cat("tau:", tau[tt], "\n")
  repeat{ 
    if (gamma_diff < eps && alpha_diff<eps | a > maxiter) {break}
    ##### Initial
    mu_i=X%*%beta
    sigma_i=exp((Z%*%gamma)/2)
    lambda_i=V%*%alpha
    r=y-mu_i
    t=r/sigma_i
    delta_i=lambda_i/sqrt(1+(lambda_i)^2)
    kappa_i2=((sigma_i)^2)*(1-(delta_i)^2)
    zeta_i=sigma_i*delta_i
    
    ### Lq weights
    WLq=diag(as.vector((2*dnorm(t)*pnorm(lambda_i*t))^(1-Lq)))
    #WLq<-diag(one) #for ML
    
    ### EM ALGORITHM WEIGHTS
    
    wsn=dnorm(lambda_i*t)/pnorm(lambda_i*t)
    u1=delta_i*t+sqrt(1-delta_i^2)*wsn   #EM U1
    u2=(1-delta_i^2)+delta_i*t*u1        #EM U2
    
    ####auxiliary matrices
    Rkappa1=diag(as.vector(1/kappa_i2))
    Rkappa2=diag(as.vector(zeta_i/kappa_i2))
    P=diag(as.vector(zeta_i*u1))
    
    R12=diag(as.vector(lambda_i/sigma_i^2))
    P12=diag(as.vector((1+2*lambda_i^2)*r/(sigma_i*sqrt(1+lambda_i^2))))
    gn=lambda_i/(1+lambda_i^2)
    one=rep(1,n)
    # LASSO ve BRIDGE
    if (penalty == 1){
      #g = 3  # LASSO iC'in 1 olacak
      diag(W_gamma) <- g*abs(gamma)^(g-2)
      W_gamma[which(!is.finite(W_gamma))] <- 0
    }
    # # SCAD   a = 3.7
    if (penalty == 2){
      diag(W_gamma) <- SCAD.derivative(abs(gamma), tau[tt])/abs(gamma)
      W_gamma[!is.finite(W_gamma)] <- 0
    }
    
    # First and second derivatives of Lq-likelihood
    U_gamma12=Reduce(`+`, lapply(1:n, function (i) WLq[i,i]^(1/2)*(-0.5+0.5*(r[i]^2/kappa_i2[i])-0.5*zeta_i[i]*r[i]*u1[i]/kappa_i2[i])*Z[i,]))   #First derivative of Lq-likelihood with sqrt(WLq) which is used for second derivative of Lq-likelihood
    U_gamma1=Reduce(`+`, lapply(1:n, function (i) WLq[i,i]*(-0.5+0.5*(r[i]^2/kappa_i2[i])-0.5*zeta_i[i]*r[i]*u1[i]/kappa_i2[i])*Z[i,]))   #First derivative of Lq-likelihood
    H_gamma1=Reduce(`+`, lapply(1:n, function (i) WLq[i,i]*(-0.5*(r[i]^2/kappa_i2[i])+0.5*((zeta_i[i]*r[i]*u1[i]/kappa_i2[i])))*(Z[i,])%*%t(Z[i,])))  #Second Derivative of log-likelihood
    H_gamma=(1-Lq)*t(t(U_gamma12))%*%U_gamma12+H_gamma1  #Second derivative of Lq-likelihood wrt gamma-gamma
    gamma_update=gamma+ginv(H_gamma-tau[[tt]]*W_gamma)%*%((tau[[tt]]*W_gamma)%*%gamma-(1/n)*U_gamma1)    #(1)
    for (gammasay in 1:d) {
      if (abs(gamma_update[gammasay]) < 0.05){gamma_update[gammasay] = 0} # 0.01 idi, 0.1 yaptD1m
    }
    gamma_diff <- norm(gamma0 - gamma_update)
    gamma=gamma_update
    gamma
    gamma_diff
    a=a+1
    
    # First and second derivatives of Lq-likelihood
    U_alpha12=Reduce(`+`, lapply(1:n, function (i) (WLq[i,i]^(1/2))*(delta_i[i]^2-lambda_i[i]*r[i]^2/sigma_i[i]^2+(sqrt(1+lambda_i[i]^2)+lambda_i[i]^2/sqrt(1+lambda_i[i]^2))*r[i]*u1[i]/sigma_i[i]-lambda_i[i]*u2[i])*V[i,]))
    U_alpha1=Reduce(`+`, lapply(1:n, function (i) WLq[i,i]*(lambda_i[i]/(1+(lambda_i[i])^2)-lambda_i[i]*r[i]^2/sigma_i[i]^2+(sqrt(1+lambda_i[i]^2)+lambda_i[i]^2/sqrt(1+lambda_i[i]^2))*r[i]*u1[i]/sigma_i[i]-lambda_i[i]*u2[i])*V[i,]))
    H_alpha1=Reduce(`+`, lapply(1:n, function (i) WLq[i,i]*((1-lambda_i[i]^2)/(1+lambda_i[i]^2)^2-(r[i]^2)/(sigma_i[i]^2)+((lambda_i[i]/sqrt(1+lambda_i[i]^2))+(lambda_i[i]^3+2*lambda_i[i])/(1+lambda_i[i]^2)^(3/2))*r[i]*u1[i]/sigma_i[i]-u2[i])*V[i,]%*%t(V[i,])))
    H_alpha=(1-Lq)*t(t(U_alpha12))%*%U_alpha12+H_alpha1  #Second derivative of Lq-likelihood wrt alpha-alpha
    alpha_update=alpha+ginv(H_alpha-tau[[tt]]*W_alpha)%*%((tau[[tt]]*W_alpha)%*%alpha-(1/n)*U_alpha1)     # (1)
    
    
    for (alphasay in 1:q) {
      if (abs(alpha_update[alphasay]) < 0.05){alpha_update[alphasay] = 0} # 0.01 idi, 0.1 yaptD1m
    }
    alpha_diff <- norm(alpha0 - alpha_update)
    alpha = alpha_update
    alpha
  }
  
  repeat{ 
    if (beta_diff < eps | c > maxiter) {break} 
    #Updates
    mu_i=X%*%beta
    sigma_i=exp((Z%*%gamma)/2)
    lambda_i=V%*%alpha
    r=y-mu_i
    t=r/sigma_i
    delta_i=lambda_i/sqrt(1+(lambda_i)^2)
    kappa_i2=((sigma_i)^2)*(1-(delta_i)^2)
    zeta_i=sigma_i*delta_i
    
    ### Update Lq weights
    WLq=diag(as.vector((2*dnorm(r/sigma_i)*pnorm(lambda_i*r/sigma_i))^(1-Lq)))
    
    ### Update EM ALGORITHM WEIGHTS
    t=r/sigma_i
    wsn=dnorm(lambda_i*t)/pnorm(lambda_i*t)
    u1=delta_i*t+sqrt(1-delta_i^2)*wsn
    u2=(1-delta_i^2)+delta_i*t*u1
    
    ####auxiliary matrices for beta_update
    Rkappa1=diag(as.vector(1/kappa_i2))
    Rkappa2=diag(as.vector(zeta_i/kappa_i2))
    
    beta_update<- ginv(t(X)%*%WLq%*%Rkappa1%*%X+tau[tt]*W_beta)%*%(t(X)%*%WLq%*%Rkappa1%*%y-t(X)%*%WLq%*%Rkappa2%*%u1)
    
    c = c+1
    for (betasay in 1:p) {
      if (abs(beta_update[betasay]) < 0.05){beta_update[betasay] = 0} 
    }
    
    beta_diff <- norm(beta0 - beta_update)
    beta = beta_update
    beta
    c=c+1
  }
  beta
  gamma
  alpha
  cat("lastbeta", beta_update, "\n")
  cat("lastgamma", gamma_update, "\n")
  cat("lastalpha", alpha_update, "\n")
  cat("beta_diff", beta_diff, "\n")
  
  l <- sum(2*dnorm((y-X%*%beta)/exp(Z%*%gamma)/2)*pnorm(V%*%alpha*(y-X%*%beta)/exp(Z%*%gamma/2)))
  
  df <- sum(beta_update != 0) + sum(gamma_update != 0) + sum(alpha_update != 0) 
  BIC <- -(1/n)*l + df*(log(n)/n) 
  BIC
  return(list(beta=beta_update, gamma=gamma_update, alpha=alpha_update, BIC=BIC, beta_diff=beta_diff,gamma_diff=gamma_diff, alpha_diff=alpha_diff, g = g))
}
#Example:
testing<-jlssm(beta=beta0, gamma=gamma0, alpha=alpha0, X, V, Z, y, Lq, penalty=1, g=2, tt=1)
testing





######################################################################################################################





