#Notes:as we have two species used in our analysis, for simplicity, this script will just use Parus monticolus as an example.
#For more Redundancy Analysis (RDA) details, please refer to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
library(vegan)    # Used to run RDA
#read genotype data
gbt_geno <- read_table("gbt.RDA.genotype",col_names = T) %>% as.data.frame()
env <- read.table("./gbt_env.txt", header=T)
pred <- env#choose factors
gbt.rda <- rda(gbt_geno ~ ., data=pred[,c(3,18,9,5,19)], scale=T)#rda
#The rda and dbrda used the same code to test significance and detect candidate SNPs
#Here just show the rda part
#Check for significance
signif.full <- anova.cca(gbt.rda, parallel=getOption("mc.cores"))
signif.axis <- anova.cca(gbt.rda, by="axis", parallel=getOption("mc.cores"))
#Identify candidate SNPs involved in local adaptation
load.rda <- scores(gbt.rda, choices=c(1:2), display="species")#The choices number depends on the number of significance axis
SNP_df <- as.data.frame(row.names(load.rda))
SNP_df$pos <- 1:nrow(SNP_df)
names(SNP_df) <- c("snp", "pos")
#3 standard deviation cutoff
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) #
ncand <- length(cand1) + length(cand2)
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")#  
cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)
#correlations of each candidate SNP with the environmental predictors
foo <- matrix(nrow=nrow(cand), ncol=5)
colnames(foo) <- c("bio3","bio18","bio9","bio5","bio19")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gbt_geno[,nam]
  foo[i,] <- apply(pred[,c(3,18,9,5,19)],2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)
length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
head(cand)
cand.pos <- left_join(cand, SNP_df, by = "snp")
write.table(cand.pos[,9], "./cand.pos.RDA.2axis", quote = F, row.names = F,col.names = F)
