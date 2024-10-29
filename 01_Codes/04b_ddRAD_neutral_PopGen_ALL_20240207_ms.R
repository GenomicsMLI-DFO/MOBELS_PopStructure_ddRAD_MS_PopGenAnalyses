# Info --------------------------------------------------------------------
#
# Overview: Population genetics analyses for eastern Arctic beluga Delphinapterus leucas using ddRAD neutral loci
# 
# Authors: Luca Montana and Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2024-02-07
#
#


# Housekeeping ------------------------------------------------------------

# Set seed
#set.seed(123)

# Clear workspace
rm(list = ls())

# Verify if you're in the right directory
getwd()
current.wd <- getwd()

# Libraries
library(readxl)
library(tidyverse)
library(here)
# library(ggpubr)
library(adegenet)
library(RColorBrewer)
library(pcadapt)
# BiocManager::install("qvalue")
library(qvalue)
library(eulerr)
# library(vegan)  # for allele frequency section
# library(hierfstat)
library(QuickPop)
library(dartR)  # FST + Structure
#install.packages("ParallelStructure", repos="http://R-Forge.R-project.org")
library(ParallelStructure)

# Functions
"%nin%" <- Negate("%in%")

na.gi.count <- function(gi){  # subsetting function
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){  # Function to create a list of loci, from a genind object
  # Create vectors for each loci
  MAF.res <- adegenet::minorAllele(gi)
  NA.res  <- na.gi.count(gi)
  
  # Filter by threshold
  MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
  cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
  
  NA.loc <- names(NA.res[NA.res <= NA.trs])
  cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
  
  # LOCI with both conditions
  LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
  LOCI.res %>% length()
  
  cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
  
  return(LOCI.res)
}

count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}

# Paths
plink_path <- "/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606" 




# Data --------------------------------------------------------------------

## Population data --------------------------------------------------------

pop.data <- read.csv(file = "./00_Data/02_Dataset/Beluga_ddRAD.csv")
pop.data %>% view()


## MAF05NA05 for genlight and genind objects ------------------------------

load("./00_Data/01_Filtering.ref/01_Bringloe/01i_UniqueFinal/populations.60102snps.638ind.H06.DP.single.final.recode.vcf.adegenet.Rdata")

pop(gl.final) <- data.frame(ID_GQ = indNames(gl.final)) %>% 
  left_join(pop.data) %>% pull(Region1)
table(pop(gl.final), useNA = "ifany")
# BEL CSB FRB JAM NEH NHB NHS NWH SEH SHS SLE SWH UNG 
#  43  27  16  24  32  39  27  68 124 112  23  14  89

pop(gi.final) <- data.frame(ID_GQ = indNames(gi.final)) %>% 
  left_join(pop.data) %>% pull(Region1)
table(pop(gi.final), useNA = "ifany")

# Subset Minor Allele Frequency

LOC.MAF05.NA05 <- filter.MAF.NA(gi.final, MAF.trs = 0.05, NA.trs = 0.05)  # 26019 loci
gl.final.MAF05NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF05.NA05]  # necessary for PCA with glPCA function (adegenet package)
gi.final.MAF05NA05 <- gi.final[loc = LOC.MAF05.NA05]  # necessary for PCA (adegenet package)
hist(adegenet::minorAllele(gi.final.MAF05NA05))

# Remove outlier loci

LOC.MAF05.NA05.neutral <- read.csv("./00_Data/03b_ddRAD_Bringloe/Loc.MAF05NA05_noOUTLIERS.csv", stringsAsFactors = F)[[1]]
gl.neutral.MAF05NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF05.NA05.neutral]  # necessary for PCA with glPCA function (adegenet package)
gi.neutral.MAF05NA05 <- gi.final[loc = LOC.MAF05.NA05.neutral]  # necessary for PCA (adegenet package)
hist(adegenet::minorAllele(gi.neutral.MAF05NA05))




# PCA ---------------------------------------------------------------------
# Using both neutral SNPs and SNPs under selection

# PCA functions

plot_pca_eig <- list(geom_bar(stat = "identity"), 
                     geom_line(),
                     scale_fill_viridis_c(),
                     labs(y = "% variance", title = NULL, x = "PC axis"),
                     theme_bw(),
                     theme(axis.text.x = element_blank(), 
                           panel.grid = element_blank(),
                           axis.ticks.length.y = unit(0.1, "in"),
                           axis.ticks.length.x = unit(0, "in"),
                           legend.position = "none"))

# Prepare GL objects for PCA

## All regions - all seasons

l.gl <- lapply(ls(pattern = "gl.neutral.MAF05NA05"), function(x) get(x))  # put genlight object in a list to perform PCA using a lapply
names(l.gl) <- ls(pattern = "gl.neutral.MAF05NA05")  # name genlight objects within list

# Prepare GI objects: useful for estimation of sample heterozygosity

## All regions - all seasons

l.gi <- lapply(ls(pattern = "gi.neutral.MAF05NA05"), function(x) get(x))  # put genlight object in a list to perform PCA using a lapply
names(l.gi) <- ls(pattern = "gi.neutral.MAF05NA05")  # name genlight objects within list


# PCA - individual level

l.pca <- lapply(l.gl, function(i){
  pca <- glPca(i, center = TRUE, scale = FALSE, parallel = TRUE, n.core = 20, nf = 1000)  # PCA
  pca
})

l.pca.mat <- lapply(l.pca, function(y){
  pca.mat <- data.frame(y$scores)  # score matrices
  print(y)
  pca.mat
})


if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01b_neutral_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01b_neutral_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01b_neutral_SNPs"))
}

# save(list = c("l.pca"), file = "./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_neutral.Rdata")
load("./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_neutral.Rdata")

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/", ".gitignore")) ){
  cat("*.Rdata", "*.RData", "!.gitignore", sep = "\n",
      file = file.path("./02_Results/01_ddRAD_Bringloe/01_PCA/", ".gitignore")) 
}

# Eig var

MAF05NA05.neutral.var.638 <- pca_var(l.pca[["gl.neutral.MAF05NA05"]], nInd(gl.final)-1) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.3,0.6,0.9,1.2)) +
  # geom_hline(yintercept = 0, col = "grey20") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
MAF05NA05.neutral.var.638

MAF05NA05.neutral.var.30 <- pca_var(l.pca[["gl.neutral.MAF05NA05"]], 30) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.3,0.6,0.9,1.2)) +
  # geom_hline(yintercept = 0, col = "grey20") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
MAF05NA05.neutral.var.30

# With inset
MAF05NA05.neutral.var <- MAF05NA05.neutral.var.30 + annotation_custom(ggplotGrob(MAF05NA05.neutral.var.638), xmin = 13, xmax = 30,
                                                                      ymin = 0.6, ymax = 1.18)
pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01b_neutral_SNPs/Var.PCA.global.inset.MAF05NA05.neutral.240209.pdf", width = 13, height = 10)
MAF05NA05.neutral.var
dev.off()



## Figures (MAF05NA05) ----------------------------------------------------

region.labels <- c("Belcher Islands","Cumberland Sound","Frobisher Bay","South-East Hudson Bay","James Bay","North-East Hudson Bay","North Hudson Strait",
                   "South Hudson Strait","Saint Lawrence estuary","Ungava Bay","Western Hudson Bay")
names(region.labels) <- c("BEL","CSB","FRB","SEH","JAM","NEH","NHS","SHS","SLE","UNG","WHB")

# All regions and seasons

neutral.pca <- l.pca[["gl.neutral.MAF05NA05"]] %>% QuickPop::pca_scoretable(naxe = 8) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>%
  droplevels() %>%  # to remove unused Pop levels
  mutate(Season = ifelse(Month %in% c(7,8), "Summer",
                         "Non-summer"))
table(neutral.pca$Region2)
neutral.pca$Region2 <- factor(neutral.pca$Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","SEH","BEL","JAM","WHB"))
table(neutral.pca$Season, useNA = "ifany")
neutral.pca$Season <- factor(neutral.pca$Season, levels = c("Non-summer","Summer"))


## Figures (MAF05NA05) ----------------------------------------------------

## All-in-one

gPCA.MAF05NA05.1v2 <- ggplot(data = neutral.pca, aes(x = score.PC1, y = score.PC2, col = Region2, shape = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("deepskyblue",0.65),alpha("springgreen4",0.65),alpha("plum2",0.65),alpha("purple3",0.65),
                                                         alpha("orchid3",0.65),alpha("palevioletred2",0.65),alpha("salmon",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("orange",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_shape_manual(values = c(6,19)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.MAF05NA05.1v2

gPCA.MAF05NA05.1v3 <- ggplot(data = neutral.pca, aes(x = score.PC1, y = score.PC3, col = Region2, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("deepskyblue",0.65),alpha("springgreen4",0.65),alpha("plum2",0.65),alpha("purple3",0.65),
                                                         alpha("orchid3",0.65),alpha("palevioletred2",0.65),alpha("salmon",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("orange",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC3 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[3] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.MAF05NA05.1v3
# ggsave(filename = "02_Results/01_ddRAD_Bringloe/01_PCA/01b_neutral_SNPs/PCA.neutral.MAF05NA05.1v3.20241017.png", plot = gPCA.MAF05NA05.1v3, 
#        width = 14, height = 8, unit = 'in')

gPCA.MAF05NA05.3v4 <- ggplot(data = neutral.pca, aes(x = score.PC3, y = score.PC4, col = Region2, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1) +  # added stroke for presentation
  scale_y_continuous(breaks = c(-10,-5,0,5,10,15,20,25)) +
  scale_color_manual(name = "Harvest region", values = c(alpha("deepskyblue",0.65),alpha("springgreen4",0.65),alpha("plum2",0.65),alpha("purple3",0.65),
                                                         alpha("orchid3",0.65),alpha("palevioletred2",0.65),alpha("salmon",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("orange",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # squares for spring-fall, 16 for summer, triangles for NA
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[["gl.neutral.MAF05NA05"]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.MAF05NA05.3v4


## Figures (MAF05NA05) - NA -----------------------------------------------

l.na.info <- lapply(l.gl, function(w){
  na.info <- data.frame(ID_GQ = indNames(w),
                        NNA = count.ind.na.gl(w))
  na.info
})

gPCA.NA.MAF05NA05.neutral.1v2 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[[1]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[[1]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_colourbar(title = "% NA", order = 1)) +
  theme(legend.position = c(0.2,0.25),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.neutral.1v2

gPCA.NA.MAF05NA05.neutral.1v3 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC3, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[[1]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC3 (", QuickPop::pca_var(l.pca[[1]])$p.eig[3] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_colourbar(title = "% NA", order = 1)) +
  theme(legend.position = c(0.2,0.25),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.neutral.1v3

gPCA.NA.MAF05NA05.neutral.3v4 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[[1]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[[1]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.neutral.3v4

gPCA.NA.MAF05NA05.neutral.5v6 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  labs(x = paste0("PC5 (", QuickPop::pca_var(l.pca[[1]])$p.eig[5] %>% round(3) * 100, "%)"),
       y = paste0("PC6 (", QuickPop::pca_var(l.pca[[1]])$p.eig[6] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.neutral.5v6

gPCA.NA.neutral.1v2.3v4.5v6 <- ggpubr::ggarrange(gPCA.NA.MAF05NA05.neutral.1v2,
                                         gPCA.NA.MAF05NA05.neutral.3v4,
                                         gPCA.NA.MAF05NA05.neutral.5v6,
                                         nrow = 3, ncol = 1,
                                         common.legend = T, align = "hv", legend = "right")
gPCA.NA.neutral.1v2.3v4.5v6




## Figures (MAF05NA05) - Heterozygosity -----------------------------------

l.ho <- lapply(l.gl, function(h){
  ho <- gl.report.heterozygosity(h, method='ind')  # Verify what gl.report.heterozigosity does
  ho
})

gPCA.Ho.MAF05NA05.neutral.1v2 <- neutral.pca %>% 
  left_join(l.ho[[1]], by = c("ID" = "ind.name")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, color = Ho)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[[1]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[[1]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_colourbar(title = "Observed\nheterozygosity", order = 1)) +
  theme(legend.position = c(0.2,0.25),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.neutral.1v2

gPCA.Ho.MAF05NA05.neutral.3v4 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[[1]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[[1]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.neutral.3v4


gPCA.Ho.MAF05NA05.neutral.5v6 <- neutral.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6, color = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  labs(x = paste0("PC5 (", QuickPop::pca_var(l.pca[[1]])$p.eig[5] %>% round(3) * 100, "%)"),
       y = paste0("PC6 (", QuickPop::pca_var(l.pca[[1]])$p.eig[6] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.neutral.5v6

gPCA.NNA.Ho.neutral.1v2.3v4.5v6 <- ggpubr::ggarrange(gPCA.NA.MAF05NA05.neutral.1v2,gPCA.Ho.MAF05NA05.neutral.1v2,
                                                     gPCA.NA.MAF05NA05.neutral.3v4,gPCA.Ho.MAF05NA05.neutral.3v4,
                                                     gPCA.NA.MAF05NA05.neutral.5v6,gPCA.Ho.MAF05NA05.neutral.5v6,
                                                     nrow = 3, ncol = 2,
                                                     common.legend = F)
pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01b_neutral_SNPs/PCA.neutral.MAF05NA05.NNA.HO.1v2.3v4.5v6.20240209.pdf", width = 15, height = 20)
gPCA.NNA.Ho.neutral.1v2.3v4.5v6
dev.off()
png(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01b_neutral_SNPs/PCA.neutral.MAF05NA05.NNA.HO.1v2.3v4.5v6.20240209.png", width = 15, height = 20, units = "in", res = 150, pointsize = 12)
gPCA.NNA.Ho.neutral.1v2.3v4.5v6
dev.off()


## Save results PCA -------------------------------------------------------

PCA.neutral.MAF05NA05 <- neutral.pca %>%
  left_join(l.ho[[1]], by = c("ID" = "ind.name")) %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>%
  mutate(Sex_qPCR = factor(Sex_qPCR, levels = c("F","M","NA")))

write.csv(PCA.neutral.MAF05NA05, file = "02_Results/01_ddRAD_Bringloe/01_PCA/01b_neutral_SNPs/PCA.neutral.MAF05NA05_n638.csv", row.names = F)




# Admixture ---------------------------------------------------------------

## vcf to plink -----------------------------------------------------------

cmd2b <- paste("--vcf", file.path("./00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers.recode.vcf"), 
               #"--recode",
               "--plink-tped",
               "--out", file.path("./00_Data/03b_ddRAD_Bringloe", "populations.24709snps.638ind.NoOutliers.recode")
)
cmd2b

A2b <- system2("vcftools", cmd2b, stdout=T, stderr=T)
A2b

## Make bed files ---------------------------------------------------------

cmd3b <- paste("--tfam", file.path("./00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers.recode.tfam"), 
               "--tped", file.path("./00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./00_Data/03b_ddRAD_Bringloe", "populations.24709snps.638ind.NoOutliers.recode")
               
)
cmd3b

A3b <- system2(file.path(plink_path, "plink"), cmd3b, stdout=T, stderr=T)
A3b

bed.file <- file.path(here::here(), "./00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)


## Admixture analysis -----------------------------------------------------

# set.seed(111)

for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("Bringloe.neutral.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.res <- data.frame(k = 1:10,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("Bringloe.neutral.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.res[i, "CV"] <- CV
  
}

CV.res$CV <- as.numeric(as.character(CV.res$CV))

CV.res %>% arrange(CV)

plot(CV.res$CV)

gg.CV <- CV.res %>% mutate(color = ifelse(k == 4, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", "Admixture.neutral.CV.pdf"), width = 4, height = 3.5)
# png(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", "Admixture.neutral.CV.png"), width = 4, height = 3.5, units = "in", res = 150, pointsize = 12)
gg.CV
dev.off()

k <- 6
Q.k2.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",2,".Q")))
Q.k3.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",3,".Q")))
Q.k4.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",4,".Q")))
Q.k5.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",5,".Q")))
Q.k6.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",6,".Q")))
Q.k7.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",7,".Q")))
Q.k8.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",8,".Q")))
Q.k9.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",9,".Q")))
Q.k10.neutral <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/", paste0("populations.24709snps.638ind.NoOutliers.recode.",10,".Q")))


Q.neutral <- bind_rows(cbind(fam$V1, Q.k10.neutral, K = 10),
                   cbind(fam$V1, Q.k9.neutral, K = 9),
                   cbind(fam$V1, Q.k8.neutral, K = 8),
                   cbind(fam$V1, Q.k7.neutral, K = 7),
                   cbind(fam$V1, Q.k6.neutral, K = 6),
                   cbind(fam$V1, Q.k5.neutral, K = 5),
                   cbind(fam$V1, Q.k4.neutral, K = 4),
                   cbind(fam$V1, Q.k3.neutral, K = 3),
                   cbind(fam$V1, Q.k2.neutral, K = 2))

head(Q.neutral)
names(Q.neutral) <- c("ID_GQ", paste0("Q", 1:10), "K")

gg.str.all <- Q.neutral %>% pivot_longer(cols =  paste0("Q", 1:10), names_to = "Group", values_to = "Q") %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_col() +
  facet_grid(K ~ Region2 , space = "free", scale = "free") +
  #scale_fill_brewer(palette = "Set1") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str.all


## Identify genetic clusters: K = 5 ---------------------------------------

Q.neutral <- Q.neutral %>%
  left_join(pop.data %>% select(ID_GQ,ID_rcpt,Region1,Location,Year,Month,Day,Community,Haplo)) %>%
  arrange(K, ID_GQ)
write.csv(Q.neutral, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/Beluga_ddRAD_meta_admixture_neutral_k2-10.csv", row.names = F)

## SLE
hist(Q.neutral$Q1[Q.neutral$K %in% 5], breaks = 20)  # SLE very differentiated (Q > 0.9)

## CSB
hist(Q.neutral$Q3[Q.neutral$K %in% 5], breaks = 20)  # CSB
table(Q.neutral$Region1[Q.neutral$Q3 > 0.7 & Q.neutral$K %in% 5])
Q.neutral %>% subset(K %in% 5 & Region1 %in% "CSB") %>% View()
Q.neutral %>% subset(K %in% 5) %>% View()  # 20 CSB, 1 WHB, 1 FRB (consistent with ResDoc), possible hybrids from CSB and WHB

## JAM
hist(Q.neutral$Q4[Q.neutral$K %in% 5], breaks = 20)  # JAM-BEL very differentiated
Q.neutral %>% subset(K %in% 5) %>% View()  # Q4 > 0.75 for clear JAM-BEL identification
Q.neutral %>% subset(K %in% 5 & Region1 %in% "JAM") %>% View()  # All JAM Q4 > 0.9 and LON Q4 = 0.83
Q.neutral %>% subset(K %in% 5 & Region1 %in% "BEL") %>% View()  # 4 SAN Q4 > 0.9 (+ 1 at 0.89). Overall gradient from 0.75 to 0.95. Difficult to distinguish from JAM

## LGR
hist(Q.neutral$Q5[Q.neutral$K %in% 5], breaks = 20)  # LGR very differentiated (and 2 hybrids Q5 = 0.36-37)
table(Q.neutral$Region1[Q.neutral$Q5 > 0.70 & Q.neutral$K %in% 5])  # 34 EHB, 11 NEH, 1 NHS, 11 SHS (0.74: 34 EHB, 10 NEH, 1 NHS, 7 SHS, no WHB; 0.6: 34 EHB, 11 NEH, 1 NHS, 12 SHS, 1 WHB)
Q.neutral %>% subset(K %in% 5) %>% View()  # A few individuals from NEH, SHS, and WHB betwee 0.6 and 0.74 ( I consider those above 0.74 as 'real' LGR individuals harvested during migrations)
Q.neutral %>% subset(K %in% 5 & Region1 %in% "EHB") %>% View()  #Q5 > 0.74 for clear LGR identification.

## HBSC
hist(Q.neutral$Q2[Q.neutral$K %in% 5], breaks = 20)  # HBSC

Q.res.neutral <- Q.neutral %>% subset(K %in% 5) %>% 
  mutate(Membership.neutral = ifelse(Q1 > 0.5, "SLE",
                                     ifelse(Q3 > 0.5, "CSB",
                                            ifelse(Q4 > 0.5, "JAM",
                                                   ifelse(Q5 > 0.5, "LGR",
                                                          "HBSC"))))) %>% 
  select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Community, Haplo, K, Q1, Q2, Q3, Q4, Q5, Membership.neutral)

write.csv(Q.res.neutral, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/Beluga_ddRAD_meta_admixture_neutral_K5.csv", row.names = F)




# FST ---------------------------------------------------------------------

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03b_neutral_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03b_neutral_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03b_neutral_SNPs"))
}

head(Q.res.neutral)
Q.res.neutral %>% group_by(Membership.neutral) %>% summarise(N = n())

# SLE vs CSB vs JAM vs LGR vs PAN

## All IDs

gl.neutral.5Q <- gl.neutral.MAF05NA05
pop(gl.neutral.5Q) <- data.frame(ID_GQ = indNames(gl.neutral.5Q)) %>% 
  left_join(Q.res.neutral) %>% pull(Membership.neutral)
table(pop(gl.neutral.5Q), useNA = "ifany")

FST.neutral.5Q <- gl.fst.pop(gl.neutral.5Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same

## 5 ID per cluster: all in common between lcWGS and ddRAD datasets but for LGR cluster (only two in common)

id.fst.5Q <- c(#"S_20_01266","S_20_01428","S_20_01513","S_20_01530","S_20_04119",  # RES not present in ddRAD dataset
  "S_20_01788","S_20_03480","S_20_01766","S_20_03664","S_20_01763",  # SLE
  "S_20_00622","S_20_00710","S_20_00711","S_20_03448","S_20_03453",  # JAM
  "S_20_00776","S_20_00661","S_20_00703","S_20_01418","S_20_03509",  # CSB
  "S_20_01865","S_20_00621","S_20_01639","S_20_00990","S_20_03339",  # LGR
  "S_20_03515","S_20_01019","S_20_00751","S_20_00948","S_20_02539")  # PAN

length(id.fst.5Q)
Q.res.neutral %>% subset(ID_GQ %in% id.fst.5Q) %>% nrow()
Q.res.neutral %>% subset(ID_GQ %in% id.fst.5Q) %>% group_by(Region1) %>% summarise(N = n())
Q.res.neutral %>% subset(ID_GQ %in% id.fst.5Q) %>% group_by(Membership.neutral) %>% summarise(N = n())  

gl.red.neutral.5Q <- gl.keep.ind(gl.neutral.5Q, id.fst.5Q, recalc = T, verbose = 5)  # keep only 35 individual of choice
table(pop(gl.red.neutral.5Q), useNA = "ifany")
FST.red.neutral.5Q <- gl.fst.pop(gl.red.neutral.5Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same


save(list = c("FST.neutral.5Q","FST.red.neutral.5Q"),
     file = file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_neutral_K5.Rdata"))
load(file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_neutral_K5.Rdata"))

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/03_FST/", ".gitignore")) ){
  cat("*.Rdata", "!.gitignore", sep = "\n",
      file = file.path("./02_Results/01_ddRAD_Bringloe/03_FST/", ".gitignore")) 
}


# Functions

table.fst <- function(fst){
  res <-  fst$Bootstraps %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")
  
  return(res)
  
}

heat.fst <- function(fst){
  res <- bind_rows(table.fst(fst),
                   table.fst(fst) %>% mutate(Pop1.int = Population2,
                                             Pop2.int = Population1,
                                             Population1 = Pop1.int,
                                             Population2 = Pop2.int) %>%
                     dplyr::select(-c(Pop1.int, Pop2.int))
  )
  
  return(res)
  
}


## Figures: SLE vs CSB vs JAM-BEL vs LGR vs HBSC --------------------------

FST.table.neutral.5Q <- FST.neutral.5Q %>% table.fst()
FST.table.red.neutral.5Q <- FST.red.neutral.5Q %>% table.fst()

FST.neutral.5Q %>% table.fst() %>%
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst),
            Max.Pvalue = round(max(`p-value`),))
# Mean         sd        Min        Max Max.Pvalue
# 0.04171029 0.02646162 0.01197309 0.08265666          0

FST.red.neutral.5Q %>% table.fst() %>%
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst),
            Max.Pvalue = round(max(`p-value`),))
# Mean         sd       Min       Max Max.Pvalue
# 0.05160079 0.03223427 0.0104226 0.1113164          0


FST.ddRAD.neutral.5Q <- FST.neutral.5Q %>% heat.fst() %>% 
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)
write.csv(FST.ddRAD.neutral.5Q, file = "./02_Results/01_ddRAD_Bringloe/03_FST/03b_neutral_SNPs/FST_ddRAD_neutral.csv", row.names = F)

FST.ddRAD.red.neutral.5Q <- FST.red.neutral.5Q %>% heat.fst() %>%
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)

gFST.ddRAD.neutral.5Q <-  FST.ddRAD.neutral.5Q %>%
  mutate(Population1 = factor(Population1, levels = c("HBSC","JAM","LGR","CSB","SLE")),
         Population2 = factor(Population2, levels = c("HBSC","JAM","LGR","CSB","SLE"))) %>% 
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) +
  geom_tile(colour=  "gray") +
  scale_fill_viridis_c()+
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(angle = 0),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 18, color = "black"))
gFST.ddRAD.neutral.5Q




# Heterozygosity and Fis --------------------------------------------------

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04b_neutral_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04b_neutral_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04b_neutral_SNPs"))
}

vcf.path <- file.path("00_Data/03b_ddRAD_Bringloe/populations.24709snps.638ind.NoOutliers.recode.vcf")

cmd4 <- paste("--vcf", file.path(current.wd, vcf.path),
              "--het")
cmd4

setwd(file.path(current.wd, "02_Results/01_ddRAD_Bringloe/04_Heterozygosity/04b_neutral_SNPs/"))  # specify wd so that .ldepth.mean file is written in the right place

A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
tail(A4)

cat(file = "Heterozygosity.ddRAD.log",
    "\n", cmd4, "\n",
    A4, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

het.ddRAD.neutral <- read.delim("./02_Results/01_ddRAD_Bringloe/04_Heterozygosity/04b_neutral_SNPs/out.het", skip=0, sep = "\t", header = T )
het.ddRAD.neutral %>% head()

pop.het.ddRAD.neutral <- Q.res.neutral %>% left_join(het.ddRAD.neutral, by = c("ID_GQ"="INDV")) %>%
  select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Community, Membership.neutral, O.HOM., E.HOM., N_SITES, 'F') %>% 
  rename(Fis = 'F') %>% 
  mutate(Marker = "ddRAD")

pop.het.lcWGS <- read.csv("./02_Results/01_ddRAD_Bringloe/04_Heterozygosity/beluga_lcWGS_9v23.csv", stringsAsFactors = F)
pop.het.lcWGS %>% head()
pop.het.lcWGS <- pop.het.lcWGS %>% rename(Fis = 'F') %>% 
  mutate(Membership = str_replace(Membership, "CBS", "CSB"),
         Membership = str_replace(Membership, "JAM-BEL", "JAM"),
         Membership = str_replace(Membership, "PAN", "HBSC")) %>%
  mutate(Marker = "lcWGS")

pop.het.ddRAD.neutral <- pop.het.ddRAD.neutral %>% mutate(ID = gsub("_rep", "", ID_GQ))
length(which(pop.het.ddRAD.neutral$ID %in% pop.het.lcWGS$ID))  # 73 samples in both dataset
nrow(pop.het.ddRAD.neutral) + nrow(pop.het.lcWGS) - length(which(pop.het.ddRAD.neutral$ID %in% pop.het.lcWGS$ID))  # 905

pop.het.all <- pop.het.ddRAD.neutral %>% select(ID_GQ, Membership.neutral, O.HOM., E.HOM., N_SITES, Fis, Marker) %>% rename("Membership" = "Membership.neutral") %>% 
  bind_rows(pop.het.lcWGS %>% select(ID_GQ=ID, Membership, O.HOM., E.HOM., N_SITES, Fis, Marker))  # ID_GQ=ID renames ID column in pop.het.lcWGS to ID_GQ, so the merging ends well

gFIS.all <- pop.het.all %>% subset(Membership %nin% "WHB") %>% 
  ggplot(aes(y = as.numeric(as.character(Fis)),
             x = factor(Membership, levels = c("RES","HBSC","JAM","LGR","CSB","SLE")),
             col = factor(Membership, levels = c("RES","HBSC","JAM","LGR","CSB","SLE")),
             shape = Marker)) +
  ggforce::geom_sina(alpha = 0.4, size = 8) +
  scale_color_manual(name = "Genetic cluster", values = c(alpha("skyblue2",0.8),alpha("royalblue1",0.8),alpha("orange",0.8),alpha("red3",0.8),
                                                          alpha("springgreen4",0.8),alpha("deepskyblue",0.8)),
                     labels = c("Resolute passage","Panmictic","James Bay","Little-Great Whale Rivers","Cumberland Sound","St. Lawrence estuary")) +
  scale_shape_manual(name = "Marker", values = c(19, 17)) +
  geom_violin(position = "dodge", col = "black", fill = "transparent", draw_quantiles = c(0.5)) +
  geom_hline(yintercept = c(0), lty = "dashed", col = "darkgray") +
  labs(y = "Fis", x = "") +
  theme_bw() +
  guides(colours = guide_legend(override.aes = list(size = 8)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = c(0.9,0.2),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        legend.box = "vertical") +
  theme(axis.text = element_text(size = 25, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        axis.title = element_text(size = 27, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
pdf(file = "02_Results/01_ddRAD_Bringloe/04_Heterozygosity/FIS.lcWGS.ddRAD.neutral.240209.pdf", width = 20, height = 10)
gFIS.all
dev.off()


# Correlation FIS ddRAD vs lcWGS + difference between means

pop.het.all.summary <- pop.het.all %>% subset(Membership %nin% c("WHB")) %>%
  group_by(Membership, Marker) %>% 
  summarise(Mean_Ho = mean(O.HOM.),
            Mean_He = mean(E.HOM.),
            Fis_mean = mean(Fis),
            Fis_median = round(median(Fis), digits = 3),
            Fis_SD = sd(Fis, na.rm = T))
pop.het.all.summary

cor.test(pop.het.all.summary$Fis_median[pop.het.all.summary$Marker %in% "ddRAD"],
         pop.het.all.summary$Fis_median[pop.het.all.summary$Marker %in% "lcWGS" & pop.het.all.summary$Membership %nin% "RES"], method = "s")  # rho = 0.8 (0.6 without SLE), p-value = 0.1333
plot(pop.het.all.summary$Fis_median[pop.het.all.summary$Marker %in% "ddRAD"],
     pop.het.all.summary$Fis_median[pop.het.all.summary$Marker %in% "lcWGS" & pop.het.all.summary$Membership %nin% "RES"])


wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 40332, p-value = 1.294e-07
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 1079, p-value = 0.01405
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 513, p-value = 0.3383
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 502, p-value = 0.0003404
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 392, p-value = 6.855e-10

# Mann-Whitney U test
## ddRAD
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 62, p-value = 1.665e-14
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 1580.5, p-value = 1.647e-08
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 12974, p-value = 0.2983
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 3334.5, p-value = 7.852e-13
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 3, p-value = 7.603e-16
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 376, p-value = 0.2852
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 2011, p-value = 1.328e-09
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 2, p-value = 1.894e-11
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 172, p-value = 6.665e-07
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 51, p-value = 2.549e-06

## lcWGS
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "RES" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 16, p-value = 0.0002432
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "RES" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 155, p-value = 0.4166
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "RES" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 122, p-value = 0.2264
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "RES" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 188, p-value = 0.8491
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "RES" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 1123, p-value = 0.4315
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 233, p-value = 1.544e-10
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 3177, p-value = 0.9012
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 2598, p-value = 0.272
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "HBSC" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 3823.5, p-value = 0.2489
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 57, p-value = 4.996e-09
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 674, p-value = 0.2576
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "JAM" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 532, p-value = 0.09068
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 13, p-value = 5.682e-09
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 248, p-value = 0.2737
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 15, p-valure = 1.185e-10