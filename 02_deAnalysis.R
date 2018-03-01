# File: 02_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 01/03/2018

## load the data
source('header.R')
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
g_did = 25
dbListFields(db, 'MetaFile')

dfFiles = dbGetQuery(db, 'select * from MetaFile where idData=25')

dfFiles

## load the 2 files
n = paste0(dfFiles$location, dfFiles$name)
dfCovariates = read.csv(n[1], header=T)
dfData = read.csv(n[2], header=T, row.names=1, sep = '\t')

dim(dfCovariates)
dim(dfData)

colnames(dfData)
str(dfCovariates)
table(as.character(dfCovariates$SampleName) %in% colnames(dfData))

## drop the outlier sample according to first pass analysis
## IL2217_2
dfCovariates = dfCovariates[!as.character(dfCovariates$SampleName) == 'IL2217_2',]
dfData = dfData[,!colnames(dfData) %in% 'IL2217_2']
i = match(as.character(dfCovariates$SampleName), colnames(dfData))

## there is a discrepency where we have some samples duplicated in the covariates table
## this is due to the those samples having multiple lanes but being merged at the count table
## after doing first pass analysis, drop these duplicated covariates
i = match(colnames(dfData), as.character(dfCovariates$SampleName))
dfCovariates = dfCovariates[i,]
identical(as.character(dfCovariates$SampleName), colnames(dfData))

## first normalise the data
mData = as.matrix(dfData)
dim(mData)

# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

# ## duplicate the missing samples to make even with covariates
# mData = mData[,i]
# dim(mData)
# mData.norm = mData.norm[,i]

#colnames(mData) = gsub(pattern = '\\.\\d', '', colnames(mData))
identical(colnames(mData), as.character(dfCovariates$SampleName))
identical(colnames(mData.norm), as.character(dfCovariates$SampleName))

# set the covariates
dfCovariates$PatientID = factor(dfCovariates$PatientID)
# There are two main clusters forming:
# Cluster 1 - patients from groups 1 and 2
# Cluster 2 - patients from groups 3 and 4
f = rep('Cluster2', times=nrow(dfCovariates))
f[dfCovariates$PatientID %in% c('1', '2')] = 'Cluster1'
dfCovariates$Clusters = factor(f)

## perform DE analysis
dfData.bk = dfData

## delete sample section after testing
mData.norm = round(mData.norm, 0)
set.seed(123)
i = sample(1:nrow(mData.norm), 100, replace = F)

dfData = data.frame(t(mData.norm[i,]))
dfData = stack(dfData)
dfData$fBatch = dfCovariates$Treatment
dfData$fAdjust1 = dfCovariates$PatientID
dfData$fAdjust2 = dfCovariates$Clusters
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)

dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1, dfData$Coef.adj2), ]

## setup the model
library(lme4)
fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
summary(fit.lme1)
ran = ranef(fit.lme1, condVar=F)

plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp3RandomEffects.stan')

library(MASS)
s = max(1, fitdistr(dfData$values, 'negative binomial')$estimate['size'])
## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# ## set initial values
ran = ranef(fit.lme1)
r1 = ran$Coef
r2 = ran$Coef.adj1
r3 = ran$Coef.adj2

initf = function(chain_id = 1) {
  list(sigmaRan1 = 0.5, sigmaRan2=0.5, sigmaRan3=2, rGroupsJitter1=r1, rGroupsJitter2=r2,
       rGroupsJitter3=r3, iSize=s)
}

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 Nclusters3=nlevels(dfData$Coef.adj2),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 NgroupMap3=as.numeric(dfData$Coef.adj2),
                 Ncol=1, 
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2, 
                    pars=c('sigmaRan1', 'sigmaRan2', 'sigmaRan3', 'betas',
                           'iSize',  
                           'rGroupsJitter1', 
                           'rGroupsJitter2',
                           'rGroupsJitter3'),
                    cores=2, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 15))
save(fit.stan, file='temp/fit.stan.norm.rds')
gc(reset = T)
