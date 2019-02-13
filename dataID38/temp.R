mFitted = extract(fit.stan)$mu
mScale = extract(fit.stan)$iSize

for (x in 1:100){
i = sample(1:nrow(mFitted), size=1)
p1 = mFitted[i, ]
p2 = mScale
s1 = rnbinom(nrow(dfData), size=p2, mu=p1)
lines(density(s1), lwd=0.5)
}
hist(s1)


## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
# The procedure for carrying out a posterior predictive model check requires specifying a test
# quantity, T (y) or T (y, θ), and an appropriate predictive distribution for the replications
# y rep [Gelman 2008]
## variance
T1_var = function(Y) return(var(Y))
## is the model adequate except for the extreme tails
T1_symmetry = function(Y, th){
  Yq = quantile(Y, c(0.90, 0.10))
  return(abs(Yq[1]-th) - abs(Yq[2]-th))
} 

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
mChecks = matrix(NA, nrow=5, ncol=3)
rownames(mChecks) = c('Variance', 'Symmetry', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('Normal', 'NormalCont', 'T')
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = nrow(dfData), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=2)
colnames(mThetas) = c('mu', 'sd')

for (i in 1:200){
  p = sample(1:nrow(mFitted), size=1)
  p1 = mFitted[p, ]
  p2 = mScale
  s1 = rnbinom(nrow(mDraws), size=p2, mu=p1)
  mDraws[,i] = s1
  #mThetas[i,] = c(p1, s)
}

### get the test quantity from the test function
t1 = apply(mDraws, 2, T1_var)
hist(t1, xlab='Test Quantity - Variance (Normal Model)', main='', breaks=50)
abline(v = var(dfData$values), lwd=2)
mChecks['Variance', 1] = getPValue(t1, var(dfData$values))
# 0.48, the result from Figure 6.4 Gelman [2013]
# The sample variance does not make a good test statistic because it is a sufficient statistic of
# the model and thus, in the absence of an informative prior distribution, the posterior
# distribution will automatically be centered near the observed value. We are not at all
# surprised to find an estimated p-value close to 1/2 . [Gelman 2008]

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, xlim=c(-12, 12), ylim=c(-12, 12), pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check (Normal Model)')
abline(0,1)
mChecks['Symmetry', 1] = getPValue(t1, t2) # we should see somewhere around 0.1 to 0.2 on repeated simulations
# The estimated p-value is 0.26, implying that any observed asymmetry in the middle of the distribution can easily be
# explained by sampling variation. [Gelman 2008]

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(dfData$values)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(dfData$values)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(dfData$values)
mChecks['Mean', 1] = getPValue(t1, t2)




# temp2$coefficients
# Estimate Std. Error    z value  Pr(>|z|)
# (Intercept)    4.83498859 0.09764307 49.5169647 0.0000000
# treatmentIL13 -0.01139021 0.10820102 -0.1052690 0.9161624
# treatmentIL9   0.01428245 0.10786401  0.1324116 0.8946587
# adjustmentSH  -0.11084446 0.10767346 -1.0294501 0.3032682
# adjustmentZB  -0.10026471 0.10752888 -0.9324445 0.3511069
# > class(temp2$coefficients)
# [1] "matrix"
# > pnorm(-abs(-0.01139021/0.10820102))*2
# [1] 0.9161624
# > 

###############################################################################
## nb glm
#### define function to optimize
library(LearnBayes)
library(MASS)
library(numDeriv)
## define log posterior
mylogpost = function(theta, data){
  size = exp(theta['size'])
  betas = theta[-1]
  if (size < 1) return(-Inf)
  ## model matrix to get fitted value
  mModMatrix = data$mModMatrix
  y = data$resp
  ## likelihood function
  lf = function(dat, pred){
    return(log(dnbinom(dat, size = size, mu = pred)))
  }
  mu = exp(mModMatrix %*% betas)
  val = sum(lf(y, mu))
  ret = val + dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 2, log=T))
  return(ret)
}

getOptimizedPValue = function(obj){
  se = sqrt(abs(diag(obj$var)))[c('beta0', 'beta1', 'beta2')]
  m = obj$mode[c('beta0', 'beta1', 'beta2')]
  pnorm(-abs(m/se))*2
}

getOptimizedSummary = function(obj){
  se = sqrt(abs(diag(obj$var)))['beta2']
  m = obj$mode['beta2']
  p = pnorm(-abs(m/se))*2
  ret  = c('Coef' = round(m, 3), 'SE'=round(se, 3), 'P-Value'=signif(p, 3))
  names(ret) = c('Coef', 'SE', 'P-Value')
  return(ret)
}

############################### fit the model on multiple genes
## calculate the distribution for sizes
# ivSizes = apply(mDat, 1, function(x){
#   getsize = function(X) fitdistr(as.integer(X), 'negative binomial')$estimate['size'];
#   tryCatch(expr = getsize(x), error=function(e) NULL)
# })
# ivSizes = unlist(ivSizes)
# 
# iSizeParam = unlist(gammaShRaFromModeSD(sd(ivSizes)/2, 2*sd(ivSizes)))
mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}


modelFunction = function(dat){
  ## prepare data before fitting
  dfData = data.frame(resp=as.integer(mDat[dat,]), pred=fCondition)
  temp = fitdistr(dfData$resp, 'negative binomial')$estimate['size']
  names(temp) = NULL
  # set starting values for optimizer and 
  # set parameters for optimizer
  lData = list(resp=dfData$resp, pred=dfData$pred)
  # hyper-hyperparameters for deflection variance, calculated using the data
  #lData$betaParam = unlist(gammaShRaFromModeSD(mode = log(sd(dfData$resp+0.5)/2), sd = log(2*sd(dfData$resp+0.5))))
  # use a jeffery's prior 
  lData$betaParam = c('shape'=0.5, 'rate'=0.0001)
  # # gamma hyperparameter for the size 
  # lData$sizeParam = iSizeParam
  start = c('size'=log(max(temp, 1)), 'beta0'=log(mean(dfData$resp)), 'beta1'=0, 'beta2' = 0, 
            'betaSigma'=log(log(sd(dfData$resp))))
  #op = optim(start, mylogpost, control = list(fnscale = -1, maxit=1000), data=lData)
  # see results of optimiser
  #if (op$convergence) cat('Optimizer convergence failed')
  ## you can see the starting values
  #start2 = op$par
  #names(start2) = names(start)
  mylaplace(mylogpost, start, lData)
}



###############################################################################








