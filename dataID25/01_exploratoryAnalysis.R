# File: 01_exploratoryAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 24/01/2018

## load the data
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
# dfCovariates = dfCovariates[!as.character(dfCovariates$SampleName) == 'IL2217_2',]
# dfData = dfData[,!colnames(dfData) %in% 'IL2217_2']
# i = match(as.character(dfCovariates$SampleName), colnames(dfData))

## there is a discrepency where we have some samples duplicated in the covariates table
## this is due to the those samples having multiple lanes but being merged at the count table
## after doing first pass analysis, drop these duplicated covariates
i = match(colnames(dfData), as.character(dfCovariates$SampleName))
dfCovariates = dfCovariates[i,]

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

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData.norm+1), 'Normalised')
oDiag.2 = CDiagnosticPlots(log(mData+1), 'Original')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
# e.g. in this case it is different lanes/machines
fBatch = factor(dfCovariates$PatientID)
fBatch = factor(dfCovariates$Treatment)
fBatch = factor(dfCovariates$Lane)
fBatch = factor(dfCovariates$RunID)

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

## extreme values
oDiag.1 = calculateExtremeValues(oDiag.1)
oDiag.2 = calculateExtremeValues(oDiag.2)
m1 = mGetExtremeValues(oDiag.1)
m2 = mGetExtremeValues(oDiag.2)

## samples with most extreme values
apply(m1, 2, function(x) sum(x > 0))
apply(m2, 2, function(x) sum(x > 0))

## variables that are contributing to this
v1 = apply(m1, 1, function(x) sum(x > 0))
v2 = apply(m2, 1, function(x) sum(x > 0))

which(v1 > 0)
which(v2 > 0)



