# File: 03_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams
# Date: 07/03/2018


lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.01,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.01,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.01,]
dfContrast4.sub = ldfData[[4]][ldfData[[4]]$adj.P.Val < 0.01,]
dfContrast5.sub = ldfData[[5]][ldfData[[5]]$adj.P.Val < 0.01,]

library(VennDiagram)

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub), rownames(dfContrast4.sub),
             rownames(dfContrast5.sub))
names(ldfData)
names(lVenn) = gsub('results//DEAnalysis(\\w+)VsControl.xls', '\\1', names(ldfData))
# calculate overlaps
#lVenn.overlap = calculate.overlap(lVenn)
venn.diagram(lVenn, filename = 'results/venn_all_contrasts.tif', margin=0.1)
