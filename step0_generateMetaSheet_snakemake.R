##################################################
## Project: macaque LMD bulk RNA-seq
## Script purpose: generate metasheet for snakemake workflow
## Date: 2020-06-23
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("~/lustrelyt/SP/bulkRNA")


## Section: load data
##################################################
df <- read.csv('data/metasheet.csv')
df <- df[,1:3]
region.ls <- unique(df$Regions)
coln.ls <- paste0('comp_',region.ls,'vsOther')
df <- data.frame(cbind(df, matrix(1, nrow=32, ncol=8)))
colnames(df)[4:11] <- coln.ls
for(i in 4:11){
  df[which(df$Regions == region.ls[i-3]), i] <- 2
}
# 1 means control   2 means treatment 
# comparison column names start with 'comp_'
write.csv(df, file = 'data/metasheet.csv', col.names = T, row.names = F, quote = F)
