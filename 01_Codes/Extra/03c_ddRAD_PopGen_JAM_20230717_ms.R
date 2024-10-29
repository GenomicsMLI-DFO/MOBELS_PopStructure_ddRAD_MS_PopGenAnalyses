# Info --------------------------------------------------------------------
#
# Overview: Population genetics analyses for eastern Arctic beluga Delphinapterus leucas using ddRAD loci - PAN cluster
# 
# Author: Luca Montana and Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-07-17
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
library(adegenet)
library(RColorBrewer)
library(qvalue)
library(eulerr)
library(QuickPop)
library(dartR)  # FST + Structure
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

# Paths
plink_path <- "/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606" 




# Data --------------------------------------------------------------------

pop.data <- read.csv(file = "./00_Data/02_Dataset/Beluga_ddRAD.csv")
adx.all <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_ALL_K5.csv", stringsAsFactors = F)
jam.ids <- adx.all %>% subset(Membership.all %in% "JAM") %>% pull(ID_GQ)
jam.ids  # N = 41

jam.data <- pop.data %>% filter(ID_GQ %in% jam.ids)



## ddRAD ------------------------------------------------------------------

if(!file.exists(file.path("./00_Data", "03_ddRAD_Bringloe"))){
  dir.create(file.path("./00_Data", "03_ddRAD_Bringloe"))
  print(file.path("./00_Data", "03_ddRAD_Bringloe"))
}

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
maf <- 0.05
LOC.MAF05.NA05 <- filter.MAF.NA(gi.final, MAF.trs = maf, NA.trs = 0.05)  # 26019 loci

gl.final.MAF05NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF05.NA05]
gi.final.MAF05NA05 <- gi.final[loc = LOC.MAF05.NA05]
hist(adegenet::minorAllele(gi.final.MAF05NA05))

# Filter for JAM individuals only
gl.jam.MAF05NA05 <- gl.keep.ind(gl.final.MAF05NA05, jam.ids, recalc = T, verbose = 5)  # keep only 35 individual of choice
table(pop(gl.jam.MAF05NA05))
# BEL JAM SEH 
#  18  19   4
gi.jam.MAF05NA05 <- gi.final.MAF05NA05[indNames(gi.final.MAF05NA05) %in% jam.ids]
table(pop(gi.jam.MAF05NA05))




# PCA ---------------------------------------------------------------------

PCA.JAM <- glPca(gl.jam.MAF05NA05, center = TRUE, scale = FALSE, parallel = TRUE, n.core = 20, nf = 1000)  # PCA
PCA.mat.JAM <- data.frame(PCA.JAM$scores)

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01d_all_SNPs_JAM"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01d_all_SNPs_JAM"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01d_all_SNPs_JAM"))
}

# save(PCA.JAM, file = "./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_jam.Rdata")
load("./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_jam.Rdata")

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/", ".gitignore")) ){
  cat("*.Rdata", "*.RData", "!.gitignore", sep = "\n",
      file = file.path("./02_Results/01_ddRAD_Bringloe/01_PCA/", ".gitignore")) 
}

# Eig var

PCA.JAM.var.41 <- pca_var(PCA.JAM, nInd(gl.jam.MAF05NA05)-1) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,5),breaks=c(0,1,2,3,4,5)) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
PCA.JAM.var.41

PCA.JAM.var.10 <- pca_var(PCA.JAM, 10) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,5),breaks=c(0,1,2,3,4,5)) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
PCA.JAM.var.10

# With inset
PCA.JAM.var <- PCA.JAM.var.10 + annotation_custom(ggplotGrob(PCA.JAM.var.41), xmin = 5, xmax = 10,
                                                  ymin = 3.2, ymax = 4.8)
PCA.JAM.var


## Figures (MAF05NA05) ----------------------------------------------------

region.labels <- c("Belcher Islands","South-East Hudson Bay","James Bay")
names(region.labels) <- c("BEL","SEH","JAM")

# All regions and seasons

JAM.pca <- PCA.JAM %>% QuickPop::pca_scoretable(naxe = 8) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>%
  droplevels() %>%  # to remove unused Pop levels
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))))


table(JAM.pca$Region2)
JAM.pca$Region2 <- factor(JAM.pca$Region2, levels = c("SEH","BEL","JAM"))
table(JAM.pca$Season, useNA = "ifany")
JAM.pca$Season <- factor(JAM.pca$Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))

gPCA.JAM.MAF05NA05.1v2 <- ggplot(data = JAM.pca, aes(x = score.PC1, y = score.PC2, col = Region2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1, size = 7) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("red3",0.65),alpha("chocolate3",0.65),alpha("orange",0.65)),
                     labels = c("South-East Hudson Bay","Belcher Islands","James Bay")) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(PCA.JAM)$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(PCA.JAM)$p.eig[2] %>% round(3) * 100, "%)"),
       tag = "B") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(legend.position = c(0.725,0.825),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal",
        legend.background = element_blank()) +  # transparent background legend
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.tag.position = c(0.15,1.05))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01d_all_SNPs_JAM/PCA.JAM.MAF05NA05.1v2.20230717.pdf", width = 13, height = 10)
gPCA.JAM.MAF05NA05.1v2
dev.off()

gPCA.JAM.MAF05NA05.3v4 <- ggplot(data = JAM.pca, aes(x = score.PC3, y = score.PC4, col = Region2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1, size = 7) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("red3",0.65),alpha("chocolate3",0.65),alpha("orange",0.65)),
                     labels = c("South-East Hudson Bay","Belcher Islands","James Bay")) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(PCA.JAM)$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(PCA.JAM)$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01d_all_SNPs_JAM/PCA.JAM.MAF05NA05.3v4.20230717.pdf", width = 13, height = 10)
gPCA.JAM.MAF05NA05.3v4
dev.off()


## Save results PCA -------------------------------------------------------

write.csv(JAM.pca, file = "./00_Data/03_ddRAD_Bringloe/01d_all_SNPs_JAM/PCA.jam.MAF05NA05_n41.csv", row.names = F)




# Admixture ---------------------------------------------------------------

jam.ids  # N = 41

# Filter by specimen ID: only PAN individuals

cmd4a <- paste("--vcf", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.vcf"),
               "--recode",
               paste("--indv", jam.ids, collapse = " "),
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", paste0("populations.26019snps.",length(jam.ids),"ind.jam"))
)
cmd4a

A4a <- system2("vcftools", cmd4a, stdout=T, stderr=T)
tail(A4a)

# Make tped tfam files

cmd4b <- paste("--vcf", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.41ind.jam.recode.vcf"), 
               #"--recode",
               "--plink-tped",
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.41ind.jam.recode")
)
cmd4b

A2b <- system2("vcftools", cmd4b, stdout=T, stderr=T)
A2b

# Make bed file

cmd4c <- paste("--tfam", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.41ind.jam.recode.tfam"), 
               "--tped", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.41ind.jam.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.41ind.jam.recode")
               
)
cmd4c

A4c <- system2(file.path(plink_path, "plink"), cmd4c, stdout=T, stderr=T)
A4c


# Admixture analysis

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD_Bringloe/populations.26019snps.41ind.jam.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)

for(k in 1:4){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD_Bringloe/02_Admixture/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("JB.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.jam <- data.frame(k = 1:4,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.jam)){
  # Which k
  k <- CV.jam[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/", paste0("JB.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.jam[i, "CV"] <- CV
  
}

CV.jam$CV <- as.numeric(as.character(CV.jam$CV))

CV.jam %>% arrange(CV)

plot(CV.jam$CV)

gg.CV.jam <- CV.jam %>% mutate(color = ifelse(k == 1, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/", "Admixture.JAM.CV.pdf"), width = 4, height = 3.5)
gg.CV.jam
dev.off()

k <- 4
Q.k2.jam <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/", paste0("populations.26019snps.41ind.jam.recode.",2,".Q")))
Q.k3.jam <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/", paste0("populations.26019snps.41ind.jam.recode.",3,".Q")))
Q.k4.jam <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/", paste0("populations.26019snps.41ind.jam.recode.",4,".Q")))

Q.jam <- bind_rows(cbind(fam$V1, Q.k4.jam, K = 4),
                   cbind(fam$V1, Q.k3.jam, K = 3),
                   cbind(fam$V1, Q.k2.jam, K = 2))

head(Q.jam)

names(Q.jam) <- c("ID_GQ", paste0("Q", 1:4), "K")
#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str.jam <- Q.jam %>% pivot_longer(cols =  paste0("Q", 1:4), names_to = "Group", values_to = "Q") %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_col() +
  facet_grid(K ~ Region1, space = "free", scale = "free") +
  # scale_fill_manual(values = c("#00008bff","#436eeeff","#aa72dfff","#f4c9f4ff")) +
  scale_fill_manual(values = c("orange","chocolate3","royalblue1","purple")) +
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
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt")
  )
gg.str.jam

# Figure: ordered by region

d.str.jam <- Q.jam %>% pivot_longer(cols =  paste0("Q", 1:4), names_to = "Group", values_to = "Q") %>%
  left_join(pop.data %>% select(ID_GQ, Region1, Location, Year, Month)) %>% 
  mutate(Region1.full = ifelse(Region1 %in% "JAM", "James Bay",
                               ifelse(Region1 %in% "BEL", "Belcher Islands",
                                      "South-East Hudson Bay")),
         Region1.full = factor(Region1.full, levels = c("James Bay","Belcher Islands","South-East Hudson Bay")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>%
  mutate(Group.order = ifelse(is.na(Q), Group,
                              ifelse(K %in% 2, Group,
                                     ifelse(K %in% 3 & Group %in% "Q1", "Q2",
                                            ifelse(K %in% 3 & Group %in% "Q2", "Q1",
                                                   ifelse(K %in% 4 & Group %in% "Q1", "Q3",
                                                          ifelse(K %in% 4 & Group %in% "Q2", "Q1",
                                                                 ifelse(K %in% 4 & Group %in% "Q3", "Q2",
                                                                        Group))))))))

gg.str.jam <- d.str.jam %>% ggplot(aes(x = ID_GQ, y = Q, fill = Group.order)) + 
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(K ~ Region1.full, space = "free", scale = "free") +
  scale_fill_manual(values = c("orange","chocolate3","moccasin","chocolate4")) +
  # scale_fill_brewer(palette = "Set1") +
  labs(y="Membership probability",
       tag = "A") +
  theme_minimal(base_size = 20, base_family = "Helvetica") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        # axis.ticks.y = element_line(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 18),
        strip.text.y = element_text(angle = 0, size = 12),
        # panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
        plot.tag.position = c(0.1,0.7))
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/", "Admixture_JAM_k2_to_k4.pdf"), width = 12, height = 6)
# png(file = "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/Admixture_JAM_k2_to_k4.png", width = 12, height = 9, units = "in", res = 150, pointsize = 12)
gg.str.jam
dev.off()

gg.str.pca.jam <- cowplot::plot_grid(NULL,
                                     gPCA.JAM.MAF05NA05.1v2,
                                     NULL,
                                     nrow = 3, ncol = 1, rel_heights = c(0.95,1.5,0.05))
gg.str.pca.admx.jam <- cowplot::plot_grid(gg.str.jam,
                                          gg.str.pca.jam,
                                          nrow = 1, ncol = 2, rel_widths = c(1.75,1))
# png(file = "02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/Admixture_PCA_JAM_k2_to_k4.png", width = 18, height = 9, units = "in", res = 150, pointsize = 12)
gg.str.pca.admx.jam
dev.off()

pop.data %>% subset(ID_GQ %in% Q.jam$ID_GQ) %>% left_join(d.str.jam) %>% subset(!duplicated(ID_GQ)) %>% group_by(Region1, Season) %>% summarise(N = n()) %>% print(n = 100)
pop.data %>% select(ID_GQ,Sex_qPCR,Region1,Location,Year,Month) %>% subset(ID_GQ %in% Q.jam$ID_GQ) %>% left_join(Q.jam) %>% subset(K %in% 2) %>%
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% arrange(Region1, Month) %>% View()
pop.data %>% subset(ID_GQ %in% Q.jam$ID_GQ) %>% left_join(Q.jam %>% filter(K %in% 2)) %>%  
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))),
         Season = factor(Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))) %>% 
  group_by(Region, Location, Season, Q2 > 0.7) %>% summarise(N = n()) %>% print(n = 100)

Q.jam.meta <- Q.jam %>% 
  left_join(pop.data %>% select(ID_GQ,Region1,Location,Year,Month,Day)) %>% 
  arrange(K, ID_GQ)
write.csv(Q.jam.meta, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/Beluga_ddRAD_meta_admixture_JAM_k2-4.csv", row.names = F)


## Identify genetic clusters ----------------------------------------------

# K = 2
Q.jam

hist(Q.jam$Q1[Q.jam$K %in% 2], breaks = 20)  # JAM cluster
table(Q.jam.meta$Region1[Q.jam.meta$Q1 > 0.5 & Q.jam.meta$K %in% 2])  # 3 BEL 15 JAM 4 SEH
table(Q.jam.meta$Region1[Q.jam.meta$Q1 > 0.5 & Q.jam.meta$K %in% 2 & Q.jam.meta$Month %in% c(7,8)])  # 14 JAM

hist(Q.jam$Q2[Q.jam$K %in% 2], breaks = 20)  # BEL cluster
table(Q.jam.meta$Region1[Q.jam.meta$Q2 > 0.5 & Q.jam.meta$K %in% 2])  # 15 BEL 4 JAM
table(Q.jam.meta$Region1[Q.jam.meta$Q2 > 0.5 & Q.jam.meta$K %in% 2 & Q.jam.meta$Month %in% c(7,8)])  # 4 JAM (likely a JAM family...)

# Clusters seems much more defined than within PAN

Q.jam %>% filter(K %in% 2) %>% nrow()  # 494 belugas
id.jam <- Q.jam %>% filter(Q1 > 0.5, K %in% 2) %>% pull(ID_GQ)  # 22
id.bel <- Q.jam %>% filter(Q2 > 0.5, K %in% 2) %>% pull(ID_GQ)  # 19
id <- c(id.jam, id.bel)  # 41
Q.jam %>% filter(K %in% 2) %>% nrow()  # 41

Q.jam.meta %>% filter(K %in% 2) %>% group_by(Region1) %>% summarise(Mean_Q1 = mean(Q1),
                                                                    Mean_Q2 = mean(Q2))
# Region1 Mean_Q1 Mean_Q2
# BEL       0.191   0.809
# JAM       0.730   0.270
# SEH       0.803   0.197

Q.res.jam <- Q.jam.meta %>% filter(K %in% 2) %>% 
  mutate(Membership.jam = ifelse(Q1 > 0.5, "JAM","BEL"),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))))

write.csv(Q.res.jam, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02d_all_SNPs_JAM/Beluga_ddRAD_meta_admixture_JAM_K2.csv", row.names = F)




# FST ---------------------------------------------------------------------
# FST within JAM cluster estimated using Regions instead of clusters, since they are of difficult definition: JAM vs BEL
# (actually not so difficult to define, but are they real?)

head(jam.data)
jam.data %>% group_by(Region1) %>% summarise(N = n()) %>% print(n = 50)
jam.data <- jam.data %>% mutate(Region.jam = ifelse(Region1 %in% c("JAM","SEH"), "JAM","BEL"))

# SEH vs WHB vs else

## All IDs

gl.jam.2Q <- gl.jam.MAF05NA05
pop(gl.jam.2Q) <- data.frame(ID_GQ = indNames(gl.jam.2Q)) %>%
  left_join(jam.data) %>% pull(Region.jam)
table(pop(gl.jam.2Q), useNA = "ifany")

FST.2Q <- gl.fst.pop(gl.jam.2Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same

# save(list = c("FST.2Q"),
#      file = file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_JAM.Rdata"))
load(file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_JAM.Rdata"))


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


## Figures: JAM vs BEL ----------------------------------------------------

FST.table.2Q <- FST.2Q %>% table.fst()
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value      Fst
#         JAM         BEL           0.01432025           0.01586006       0 0.015059

FST.ddRAD.2Q <- FST.2Q %>% heat.fst() %>%
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)
