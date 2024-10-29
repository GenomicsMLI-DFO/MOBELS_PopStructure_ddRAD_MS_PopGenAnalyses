# Info --------------------------------------------------------------------
#
# Overview: Outlier loci supplementary analysis
# 
# Authors: Luca Montana and Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-08-01
#
# Overview: Detect outlier loci with pcadapt (Privé et al. 2020)



# Library -----------------------------------------------------------------

library(tidyverse)
library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")



# PCadapt -----------------------------------------------------------------

# Read .bed in PCAadapt
pcadapt.final.genotype  <- read.pcadapt("00_Data/03b_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.bed",
                                        type = "bed")  # MAF05NA05 created in 03_ddRaD_PopGen_ALL_20230614.R


pcadapt.final.snp <- read.delim("00_Data/03b_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.bim",
                                header = F) %>% pull(V2)  # MAF05NA05 created in 03_ddRaD_PopGen_ALL_20230614.R

# Run pcadapt

K.init <- 10

pcadapt.final <- pcadapt(pcadapt.final.genotype, K =K.init, min.maf = 0.01)

pcadapt.final$maf %>% min()

# Check screeplot

plot(pcadapt.final, option = "screeplot") 

# Check structure

plot(pcadapt.final, option = "scores") 

# K = 4 based on the Scree plot

pcadapt.final.k4 <- pcadapt(pcadapt.final.genotype , K = 4,  min.maf = 0.01)

plot(pcadapt.final.k4, option = "manhattan")
hist(pcadapt.final.k4$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(pcadapt.final.k4, option = "qqplot")

# Statistics
#x$pvalues 
alpha <- 0.05


qval.final.k4 <- qvalue::qvalue(pcadapt.final.k4$pvalues)$qvalues
outliers.final.k4 <-  which(qval.final.k4 < alpha)#, which(is.na(qval.final.k3))) %>% unique()

length(outliers.final.k4)

pcadapt.outliers <- pcadapt.final.snp[outliers.final.k4]


# VCF without them --------------------------------------------------------

library(Hmisc)

snp.info <- read.csv("00_Data/03b_ddRAD_Bringloe/Loc.MAF05NA05.csv")

snp.info %>% nrow()


snp.new <- snp.info %>% dplyr::filter(ID %nin% pcadapt.outliers )

#Check that it fits
nrow(snp.new) + length(pcadapt.outliers) == snp.info %>% nrow()

write.csv(snp.new  %>% select(ID), "00_Data/03b_ddRAD_Bringloe/Loc.MAF05NA05_noOUTLIERS.csv", 
          row.names = F, quote = F)


cmd <- paste("--vcf", "00_Data/03b_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.vcf", 
               "--recode",
               "--snps", file.path("00_Data/03b_ddRAD_Bringloe/Loc.MAF05NA05_noOUTLIERS.csv"),
               "--out", "00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers"
)
cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)  
A

# We identified outlier loci using PCAdapt (v4.3.3, Privé et al. 2020), an individual-based genome-scan approach which identifies SNPs 
# significantly associated with genetic structure underlying principal components (PCs) using Mahalanobis distance. The number of PCs 
# used was K=4 based on a visual inspection of the scree plot. SNPs with a q-value < 0.05 were identified as outliers (N = 1310) and 
# removed to kept only putative neutral SNPs (N = 24,709).
