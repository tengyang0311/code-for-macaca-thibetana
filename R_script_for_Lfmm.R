#Notes:as we have two species used in our analysis, for simplicity, this script will just use Parus monticolus as an example.
#This script used for running LFMM analysis
#Install the dependencies package
# install.packages("devtools")
#devtools::install_github("bcm-uga/LEA")
rm(list=ls())
library(LEA)#version 3.4.0
##read genotype data
gbt.geno <- read_delim("./gbt.genotype", 
	"\t", col_names = FALSE, escape_double = FALSE) %>% as.matrix()
#write.geno(gbt.geno, "gbt_genotypes.geno")#Save .geno file if needed
#Load bio factor
bio <-c("bio3","bio18","bio9","bio5","bio19")
for(i in bio){
gbt <- read.table(paste(i,".txt",sep=""),header = F) %>% as.matrix() %>% scale()
#Write.env(gbt, paste("gbt_",i,".env"))#Save .env file if needed
#Apply lfmm
mod<- lfmm2(input = gbt.geno, env = gbt, K = 3)
pv<- lfmm2.test(object = mod,
                 input = gbt.geno,
                 env = gbt,
                 linear = TRUE)
#Convert p-value to correction FDR
library(qvalue)
qval<- qvalue(pv$pvalues)$qvalues
alpha <- 0.05
outliers<- which(qval< alpha)
write.table(outliers,paste("gbt_",i,".txt",sep=""),row.names = FALSE)
}

