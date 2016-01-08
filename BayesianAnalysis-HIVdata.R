library(LearnBayes)
library("MCMCpack")
library(coda)

#read the data
hiv<-read.csv("Desktop/Raquel/R code/datacd4.csv")


hiv$x2<-factor(hiv$x2)
summary(hiv)
#sort(hiv$CD4PCT)
#sort(baseage)
#complete.cases(hiv)

#plots
hist(hiv$CD4PCT)
plot (density (hiv$CD4PCT) )
plot (sort (hiv$CD4PCT), pch=".")
plot (hiv$baseage, hiv$CD4PCT, xlab="baseage", ylab="percentage of CD4 cells")
plot (hiv$x2, hiv$CD4PCT, xlab="Recoded treatment", ylab="percentage of CD4 cells")
pairs(hiv)

#design matrix
#creating x2
#x1<-rep(1, length(hiv$CD4PCT))
#x<-cbind(x1,  hiv$baseage, hiv$x2)
#y<-(sqrt(hiv$CD4PCT))
#xtxi<-solve(t(x) %*% x)
#xtxi %*% t(x) %*% y

######################## Model I - Using non-informative prior ###############
#simple regression on CD4PCT cells using baseage and treatment as predictors
fit <- lm(sqrt(hiv$CD4PCT)~hiv$baseage+hiv$x2, data=hiv, x=TRUE, y=TRUE)
summary(fit)
#ls.out <- lsfit (fit$x, fit$y )

S=sum(fit$residual^2)
shape=fit$df.residual/2; rate=S/2
sigma2=1/rgamma(1,shape,rate)
shape1=(436-3)/2
MSE = sum(fit$residuals^2)/fit$df.residual
vbeta=vcov(fit)/MSE
#beta=rmvnorm(1,mean=fit$coef,sigma=xtxi*sigma2)

#using functions from LearnBayes package to sample betas
#in model 1, non-informative prior is assumed
theta.sample=blinreg(fit$y,fit$x,1000)

#posterior inferences in Model 1
hist(theta.sample$beta[,2])
apply(theta.sample$beta,2,quantile,c(.025,.5,.975))
plot(theta.sample$beta[,2],theta.sample$beta[,3], main="Draws of the parameter beta and gamma", xlab="beta", ylab="gamma")
plot(density(theta.sample$beta[,2]))
#posterior summaries
par(mfrow=c(2,2))
mu.sig.gibbs <- cbind(theta.sample$beta[,1],theta.sample$beta[,2], theta.sample$beta[,3], theta.sample$sigma^2)

# ... cbind() makes a matrix out of column vectors
# Set meaningful names to the columns of matrix mu.sig.gibbs:
colnames(mu.sig.gibbs) <- c('alpha', 'beta', 'gamma', 'sigma2')

# Let's reject first half of the simulations
#mu.sig.gibbs <- mu.sig.gibbs[-(1 : round(N.iter/2) ), ]
# Posterior means:
apply(mu.sig.gibbs, 2, mean)

# This is a fancy way of calculating
# mean(mu.sig.gibbs[,1]) and mean(mu.sig.gibbs[,2])
# with a single call.

# Posterior covariance:
cov(mu.sig.gibbs)

# Quantiles of the marginal posterior distributions:
apply(mu.sig.gibbs, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# You can obtain the same information with two calls:
# quantile(mu.sig.gibbs[,1], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
# quantile(mu.sig.gibbs[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
gibbs.mcmc <- mcmc(mu.sig.gibbs)
par(mfrow=c(1,1))
plot(gibbs.mcmc)
autocorr.plot(gibbs.mcmc)
effectiveSize(gibbs.mcmc)

#model checking 
yrep=matrix(0, nrow=436, ncol=1000)
mean=fit$x%*%t(theta.sample$beta)
for(i in 1:length(hiv$CD4PCT))
{
  for(j in 1:1000)
  {
    yrep[i,j]=rnorm(1, mean[i,j],theta.sample$sigma[j])
  }
}

#get the predicted draws
pred.draws=blinregpred(fit$x,theta.sample)
pred.sum=apply(pred.draws^2,2,quantile,c(.05,.95))
par(mfrow=c(1,1))
ind=1:length(hiv$CD4PCT)
matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,xlab="number of observations",ylab="Percentage of CD4PCT cells")
points(ind,hiv$CD4PCT,pch=19)
#out1=(hiv$CD4PCT>pred.sum[2,])
#out2=(hiv$CD4PCT<pred.sum[2,])
#text(ind[out1], hiv$CD4PCT[out1], label=hiv$newpid[out1], pos = 4)
#text(ind[out2], hiv$CD4PCT[out2], label=hiv$newpid[out2], pos = 4)
#posterior predictive p-value
mean.p<-rep(0, 436)
mean.p1<-rep(0, 436)
for(i in 1:length(hiv$CD4PCT))
{
  mean.p[i]= mean(c(yrep[,i]^2)>hiv$CD4PCT[i])
}
hist(mean.p ,breaks=50, main="posterior predictive p-values", xlab="p-values")
for(i in 1:length(hiv$CD4PCT))
{
  mean.p1[i]= mean(c(pred.draws[,i]^2)>hiv$CD4PCT[i])
}
hist(mean.p1, breaks=50)
plot.new()

#residual analysis
plot(sqrt(hiv$CD4PCT), fit$residual)
#plot(fit)
Residuals<-residuals(fit)
fitted<-fitted(fit)

#plot(fitted, Residuals)
plot(fitted, Residuals, xlab = "Fitted values",ylab = "Residuals", ylim = max(abs(Residuals)) * c(-1, 1))
abline(h = 0, lty = 5)

#bayesian residuals
beta=rmvnorm(1,mean=fit$coef,sigma=xtxi*sigma2)
residual<-rep(0,436)
predicted<-rep(0,436)
for (i in 1:length(hiv$CD4PCT)) {
  beta=rmvnorm(1,mean=fit$coef,sigma=xtxi*sigma2)
  predicted[i] <- beta[,1]+ beta[,2]*hiv$baseage[i]+beta[,3]*as.numeric(hiv$x2[i])      # Predicted values
  residual[i] <- sqrt(hiv$CD4PCT[i])-predicted[i]
}
#text(clouds_fitted, clouds_resid, labels = rownames(hiv$ne))

#plots of predicted values vs residuals
plot(predicted, residual)
abline(h = 0, lty = 5)


################################################### Model 2 using G-prior ##############

c.prior=c(0.1,0.5,5,2)
fit1=vector("list",4)
for (j in 1:4)
{
  prior=list(b0=c(0,0, 0), c0=c.prior[j])
  fit1[[j]]=blinreg(fit$y, fit$x, 1000, prior)
}
BETA=NULL

for (j in 1:4)
{
  s=data.frame(Prior=paste("g =",as.character(c.prior[j])), alpha=fit1[[j]]$beta[,1],beta=fit1[[j]]$beta[,2], gamma=fit1[[j]]$beta[,3])
  BETA=rbind(BETA,s)
}


library(lattice)
with(BETA,xyplot(alpha~beta|Prior,type=c("p","g")))

#posterior inference for Model 2
apply(fit1[[3]]$beta,2,quantile,c(.025,.5,.975))
#with g=5
mu.sig.gibbs1 <- cbind(fit1[[3]]$beta[,1], fit1[[3]]$beta[,2], fit1[[3]]$beta[,3],fit1[[3]]$sigma )
colnames(mu.sig.gibbs1) <- c('alpha', 'beta', 'gamma', 'sigma2')

# Let's reject first half of the simulations
# mu.sig.gibbs <- mu.sig.gibbs[-(1 : round(N.iter/2) ), ]
# Posterior means:
apply(mu.sig.gibbs1, 2, mean)

# This is a fancy way of calculating
# mean(mu.sig.gibbs[,1]) and mean(mu.sig.gibbs[,2])
# with a single call.
# Posterior covariance:
cov(mu.sig.gibbs1)

# Quantiles of the marginal posterior distributions:
apply(mu.sig.gibbs1, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# You can obtain the same information with two calls:
# quantile(mu.sig.gibbs[,1], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
# quantile(mu.sig.gibbs[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
gibbs.mcmc1 <- mcmc(mu.sig.gibbs1)
par(mfrow=c(1,1))
plot(gibbs.mcmc1)
autocorr.plot(gibbs.mcmc1)

#with g=5
pred.draws=blinregpred(fit$x,fit1[[3]])
pred.sum=apply(pred.draws^2,2,quantile,c(.05,.95))
par(mfrow=c(1,1))
ind=1:length(hiv$CD4PCT)
matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,xlab="Number of observations",ylab="Percentage of CD4 cells")
points(ind,hiv$CD4PCT,pch=19)
for(i in 1:length(hiv$CD4PCT))
{
  mean.p1[i]= mean(c(pred.draws[i,]^2)>hiv$CD4PCT[i])
}
hist(mean.p1, breaks=50, xlab="Posterior predictive p-values", main="Histogram")

for (i in 1:length(hiv$CD4PCT)) {
  fit1[[3]]=blinreg(fit$y, fit$x, 1, prior)
  predicted[i] <- fit1[[3]]$beta[,1]+ fit1[[3]]$beta[,2]*hiv$baseage[i]+fit1[[3]]$beta[,3]*as.numeric(hiv$x2[i])           # Predicted values
  residual[i] <- sqrt(hiv$CD4PCT[i])-predicted[i]
}
#text(clouds_fitted, clouds_resid, labels = rownames(hiv$ne))

#plot between predicted values vs residuals
plot(predicted, residual, xlab="fitted")
abline(h = 0, lty = 5)


######################### Model 3 - Using hierarchical bayesian model

for(i in 37:104)
{
  for(j in 1:436)
  {
    if(hiv$newpid[j]==i)
      hiv$newpid[j]=i-1
  }
}
for(i in 80:104)
{
  for(j in 1:436)
  {
    if(hiv$newpid[j]==i)
      hiv$newpid[j]=i-1
  }
}

for(i in 101:104)
{
  for(j in 1:436)
  {
    if(hiv$newpid[j]==i)
      hiv$newpid[j]=i-1
  }
}

m<-matrix(0, nrow=436, ncol=101)
p<-matrix(0, nrow=436, ncol=101)
k<-matrix(0, nrow=436, ncol=101)
for(i in 1:length(hiv$baseage))
{
  m[i,hiv$newpid[i]]<-1
}
for(i in 1:length(hiv$baseage))
{
  p[i,hiv$newpid[i]]<-hiv$baseage[i]
}

for(i in 1:length(hiv$baseage))
{
  k[i,hiv$newpid[i]]<-hiv$x2[i]
}

x<-cbind(m,p,k)

x.tran<-t(x) %*% x
tau.a<-runif(1, min = 0, max = 200)
#sigma#initial 121, 99, 150
mu1<-rep(0.1,101);mu2<-rep(0.4,101);mu3<-rep(0.6,101)
mu<-c(mu1,mu2,mu3)
sigma<-diag(mu)
sig1<-solve(sigma)
#mean
b1<-rep(3.8,101);b2<-rep(-0.12,101);b3<-rep(0.4,101)
b<-c(b1,b2,b3)
b.new<-as.matrix(b)
pvar=solve(sig1+x.tran/(2.3))
y<-as.matrix(sqrt(hiv$CD4PCT))
pmean=pvar%*%(sig1%*%b.new+t(x)%*%y/(2.3))
#beta.new<-rnorm(1, pmean, pvar)
beta.new <- rmvnorm(n=1000, mean=pmean, sigma=pvar)
beta.new1<-matrix(0, nrow=303, ncol=1000)
for(i in 1:303)
{
  for(j in 1:1000)
  {
    beta.new1[i,j]=rnorm(1, pmean[i,1], sqrt(pvar[i,i]))
  }
}


apply(beta.new, 2, quantile, probs = c(0.025,  0.5, 0.975))
plot(beta.new[,102], beta.new[,203], xlab="beta1 for patient 1", ylab="gamma1 for Patient 1")
#mvrnorm(n = 1, mu=pmean, Sigma=pvar)
mean1=x%*%t(beta.new)
#yrep<-rmvnorm(n=500, mean=x%*%t(beta.new), 2.4)
yrep1=matrix(0, nrow=436, ncol=1000)
for(i in 1:length(hiv$CD4PCT))
{
  for(j in 1:1000)
  {
    yrep1[i,j]=rnorm(1, mean1[i,j], sqrt(2.3))
  }
}

for(i in 1:length(hiv$CD4PCT))
{
  mean.p1[i]= mean(c(yrep1[i,]^2)>hiv$CD4PCT[i])
}
hist(mean.p1, breaks=50, xlab="p-values", main="posterior predictive p-values")
plot.new()

#residual analysis
plot(sqrt(hiv$CD4PCT), fit$residual)
#plot(fit)
Residuals<-residuals(fit)
fitted<-fitted(fit)

#plot(fitted, Residuals)
plot(fitted, Residuals, xlab = "Fitted values",ylab = "Residuals", ylim = max(abs(Residuals)) * c(-1, 1))
abline(h = 0, lty = 5)

#bayesian residuals
beta=rmvnorm(1,mean=fit$coef,sigma=xtxi*sigma2)
residual1<-rep(0,436)
predicted1<-rep(0,436)

for(i in 1:101)
{
  while(length(which(hiv$newpid==i))) {
    beta.new <- rmvnorm(n=1, mean=pmean, sigma=pvar)
    predicted1[m] <- beta.new[i]+ beta.new[101+i]*hiv$baseage[m]+beta.new[202+i]*as.numeric(hiv$x2[m])           # Predicted values
    residual1[m] <- sqrt(hiv$CD4PCT[m])-predicted[m]
    m=m+1
    if(m==437) break
  }
  
}

#text(clouds_fitted, clouds_resid, labels = rownames(hiv$ne))
plot(predicted1[1:436], residual1[1:436], xlab="fitted", ylab="residuals")
abline(h = 0, lty = 5)

pred.sum1=apply(yrep1^2,1,quantile,c(.05,.95))
par(mfrow=c(1,1))
ind=1:length(hiv$CD4PCT)
matplot(rbind(ind,ind),pred.sum1,type="l",lty=1,col=1,xlab="Number of observations",ylab="Percentage of CD4 cells")
points(ind,hiv$CD4PCT,pch=19)

mu.sig.gibbs3 <- cbind(beta.new[,4],beta.new[,105], beta.new[,206])
colnames(mu.sig.gibbs3) <- c('alpha3', 'beta1', 'gamma1')

# Let's reject first half of the simulations
#mu.sig.gibbs <- mu.sig.gibbs[-(1 : round(N.iter/2) ), ]
# Posterior means:
apply(mu.sig.gibbs3, 2, mean)

# This is a fancy way of calculating
# mean(mu.sig.gibbs[,1]) and mean(mu.sig.gibbs[,2])
# with a single call.
# Posterior covariance:
cov(mu.sig.gibbs3)
# Quantiles of the marginal posterior distributions:
apply(mu.sig.gibbs3, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# You can obtain the same information with two calls:
# quantile(mu.sig.gibbs[,1], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
# quantile(mu.sig.gibbs[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
gibbs.mcmc3 <- mcmc(mu.sig.gibbs3)
par(mfrow=c(1,1))
plot(gibbs.mcmc3)
autocorr.plot(gibbs.mcmc3)
for (j in 1:4)
{
  s=data.frame(Prior=paste("g =",as.character(c.prior[j])), alpha=fit1[[j]]$beta[,1],beta=fit1[[j]]$beta[,2], gamma=fit1[[j]]$beta[,3])
  BETA=rbind(BETA,s)
}

with(BETA,xyplot(alpha~beta|Prior,type=c("p","g")))
#effectiveSize(gibbs.mcmc)

