# File: 02_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 01/03/2018

## load the data
source('header.R')
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
g_did = 38
dbListFields(db, 'MetaFile')

dfFiles = dbGetQuery(db, 'select * from MetaFile where idData=38')

dfFiles

## load the 2 files
n = paste0(dfFiles$location, dfFiles$name)
dfCovariates = read.csv(n[1], header=T)
dfData = read.csv(n[2], header=T, row.names=1, sep = ',')

dim(dfCovariates)
dim(dfData)

colnames(dfData)
str(dfCovariates)
table(as.character(dfCovariates$SampleID) %in% colnames(dfData))

i = match(colnames(dfData), as.character(dfCovariates$SampleID))
dfCovariates = dfCovariates[i,]
identical(as.character(dfCovariates$SampleID), colnames(dfData))

## first normalise the data
mData = as.matrix(dfData)
dim(mData)

# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData) | inData < 1 ] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

table(ivProb < 0.4)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData), as.character(dfCovariates$SampleID))
identical(colnames(mData.norm), as.character(dfCovariates$SampleID))
## perform DE analysis
dfData.bk = dfData

## delete sample section after testing
mData.norm = round(mData.norm, 0)

set.seed(123)
i = sample(1:nrow(mData.norm), 30, replace = F)
dfData = data.frame(t(mData.norm[i,]))

#dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
dfData$fBatch = dfCovariates$Group
dfData$fAdjust1 = dfCovariates$Donor

dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)

dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1), ]
str(dfData)

# # setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef) + (1 | Coef.adj1), data=dfData)
# summary(fit.lme1)
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
#l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# ## set initial values
#ran = ranef(fit.lme1)
r1 = rep(0, nlevels(dfData$Coef))
#r2 = rep(0, nlevels(dfData$Coef.adj1))
#r3 = rep(0.01, nlevels(dfData$ind))

initf = function(chain_id = 1) {
  list(sigmaRan1 = 0.1, rGroupsJitter1=r1)
}

## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
#d = dfData[!duplicated(dfData$Coef), ]

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 #Nclusters2=nlevels(dfData$Coef.adj1),
                 #NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 Nsizes=nlevels(dfData$ind),
                 NgroupMap1=as.numeric(dfData$Coef),
                 #NgroupMap2=as.numeric(dfData$Coef.adj1),
                 #NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 NsizeMap=as.numeric(dfData$ind),
                 y=dfData$values) 
                 #gammaShape=l$shape, gammaRate=l$rate)



fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4,
                    pars=c('sigmaRan1', #'sigmaRan2',
                           'iSize', #'mu',
                           'rGroupsJitter1'),
                    cores=4, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 11))
save(fit.stan, file='temp/fit.stan.nb_test1.rds')

print(fit.stan, c('sigmaRan1', 'iSize'), digits=3)
print(fit.stan, c('rGroupsJitter1'))
traceplot(fit.stan, 'betas')
traceplot(fit.stan, c('sigmaRan2'))
traceplot(fit.stan, c('sigmaRan1'))


## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
iIntercept = as.numeric(extract(fit.stan)$betas)
## add the intercept to each random effect variable, to get the full coefficient
mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='control', deflection='IL9') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
library(org.Hs.eg.db)
df = select(org.Hs.eg.db, keys = as.character(rownames(dfResults)), columns = 'SYMBOL', keytype = 'ENSEMBL')
i = match(rownames(dfResults), df$ENSEMBL)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(rownames(dfResults), df$ENSEMBL)
## produce the plots 
f_plotVolcano(dfResults, 'TNFa vs Control', fc.lim=c(-2.5, 2.5))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'TNFa vs Control')
table(dfResults$adj.P.Val < 0.01)
## save the results 
write.csv(dfResults, file='results/DEAnalysisTNFaVsControl.xls')



######### do a comparison with deseq2
dfDesign = data.frame(Treatment = dfCovariates$Treatment, fAdjust1 = dfCovariates$PatientID, fAdjust2 = dfCovariates$Clusters,
                      row.names=colnames(mData))

oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment + fAdjust1)
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)
oRes = results(oDseq, contrast=c('Treatment', 'TNFa', 'Control'))
temp = as.data.frame(oRes)
i = match(rownames(dfResults), rownames(temp))
temp = temp[i,]
identical(rownames(dfResults), rownames(temp))
plot(dfResults$logFC, temp$log2FoldChange, pch=20)
table(oRes$padj < 0.01)