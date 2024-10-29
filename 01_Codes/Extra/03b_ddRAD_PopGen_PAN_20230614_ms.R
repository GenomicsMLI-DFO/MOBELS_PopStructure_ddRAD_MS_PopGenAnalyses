# Info --------------------------------------------------------------------
#
# Overview: Population genetics analyses for eastern Arctic beluga Delphinapterus leucas using ddRAD loci - PAN cluster
# 
# Authors: Luca Montana and Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-06-14
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
adx.all <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/Beluga_ddRAD_meta_admixture_ALL_K5.csv", stringsAsFactors = F)
pan.ids <- adx.all %>% subset(Membership.all %in% "PAN") %>% pull(ID_GQ)
pan.ids  # N = 497

pan.data <- pop.data %>% filter(ID_GQ %in% pan.ids)



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

# Filter for PAN individuals only
gl.pan.MAF05NA05 <- gl.keep.ind(gl.final.MAF05NA05, pan.ids, recalc = T, verbose = 5)  # keep only 35 individual of choice
table(pop(gl.pan.MAF05NA05))
# BEL CSB FRB JAM NEH NHB NHS NWH SEH SHS SLE SWH UNG 
#  25   7  15   5  21  39  26  67  86 101   2  14  89
gi.pan.MAF05NA05 <- gi.final.MAF05NA05[indNames(gi.final.MAF05NA05) %in% pan.ids]
table(pop(gi.pan.MAF05NA05))




# PCA ---------------------------------------------------------------------

PCA.PAN <- glPca(gl.pan.MAF05NA05, center = TRUE, scale = FALSE, parallel = TRUE, n.core = 20, nf = 1000)  # PCA
PCA.mat.PAN <- data.frame(PCA.PAN$scores)

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "01_PCA"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01c_all_SNPs_HBSC"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01c_all_SNPs_HBSC"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01c_all_SNPs_HBSC"))
}

# save(PCA.PAN, file = "./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_pan.Rdata")
load("./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_pan.Rdata")

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/", ".gitignore")) ){
  cat("*.Rdata", "*.RData", "!.gitignore", sep = "\n",
      file = file.path("./02_Results/01_ddRAD_Bringloe/01_PCA/", ".gitignore")) 
}

# Eig var

PCA.PAN.var.497 <- pca_var(PCA.PAN, nInd(gl.pan.MAF05NA05)-1) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,0.5),breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
PCA.PAN.var.497

PCA.PAN.var.30 <- pca_var(PCA.PAN, 30) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,0.5),breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
PCA.PAN.var.30

# With inset
PCA.PAN.var <- PCA.PAN.var.30 + annotation_custom(ggplotGrob(PCA.PAN.var.497), xmin = 13, xmax = 30,
                                                  ymin = 0.32, ymax = 0.5)
PCA.PAN.var


## Figures (MAF05NA05) ----------------------------------------------------

region.labels <- c("Belcher Islands","Cumberland Sound","Frobisher Bay","South-East Hudson Bay","James Bay","North-East Hudson Bay","North Hudson Strait",
                   "South Hudson Strait","Saint Lawrence estuary","Ungava Bay","Western Hudson Bay")
names(region.labels) <- c("BEL","CSB","FRB","SEH","JAM","NEH","NHS","SHS","SLE","UNG","WHB")

# All regions and seasons

PAN.pca <- PCA.PAN %>% QuickPop::pca_scoretable(naxe = 8) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>%
  droplevels() %>%  # to remove unused Pop levels
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))))


table(PAN.pca$Region2)
PAN.pca$Region2 <- factor(PAN.pca$Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","SEH","BEL","JAM","WHB"))
table(PAN.pca$Season, useNA = "ifany")
PAN.pca$Season <- factor(PAN.pca$Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))

gPCA.PAN.MAF05NA05.1v2 <- ggplot(data = PAN.pca, aes(x = score.PC1, y = score.PC2, col = Region2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1, size = 7) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),
                                                         alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("gray",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(PCA.PAN)$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(PCA.PAN)$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01c_all_SNPs_HBSC/PCA.PAN.MAF05NA05.1v2.20230717.pdf", width = 13, height = 10)
gPCA.PAN.MAF05NA05.1v2
dev.off()

gPCA.PAN.MAF05NA05.3v4 <- ggplot(data = PAN.pca, aes(x = score.PC3, y = score.PC4, col = Region2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1, size = 7) +  # added stroke for presentation
  scale_color_manual(name = "Harvest region", values = c(alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),
                                                         alpha("gray",0.65),alpha("gray",0.65),alpha("gray",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("gray",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  labs(x = paste0("PC3 (", QuickPop::pca_var(PCA.PAN)$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(PCA.PAN)$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01c_all_SNPs_HBSC/PCA.PAN.MAF05NA05.3v4.20230717.pdf", width = 13, height = 10)
gPCA.PAN.MAF05NA05.3v4
dev.off()


## Save results PCA -------------------------------------------------------

write.csv(PAN.pca, file = "./02_Results/01_ddRAD_Bringloe/01_PCA/01c_all_SNPs_HBSC/PCA.pan.MAF05NA05_n497.csv", row.names = F)




# Admixture ---------------------------------------------------------------

pan.ids  # N = 497

# Filter by specimen ID: only PAN individuals

cmd4a <- paste("--vcf", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.vcf"),
               "--recode",
               paste("--indv", pan.ids, collapse = " "),
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", paste0("populations.26019snps.",length(pan.ids),"ind.pan"))
)
cmd4a

A4a <- system2("vcftools", cmd4a, stdout=T, stderr=T)
tail(A4a)

# Make tped tfam files

cmd4b <- paste("--vcf", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.497ind.pan.recode.vcf"), 
               #"--recode",
               "--plink-tped",
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.497ind.pan.recode")
)
cmd4b

A2b <- system2("vcftools", cmd4b, stdout=T, stderr=T)
A2b

# Make bed file

cmd4c <- paste("--tfam", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.497ind.pan.recode.tfam"), 
               "--tped", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.497ind.pan.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.497ind.pan.recode")
               
)
cmd4c

A4c <- system2(file.path(plink_path, "plink"), cmd4c, stdout=T, stderr=T)
A4c


# Admixture analysis

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD_Bringloe/populations.26019snps.497ind.pan.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)

for(k in 1:4){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("Panmictic.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.pan <- data.frame(k = 1:4,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.pan)){
  # Which k
  k <- CV.pan[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC", paste0("Panmictic.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.pan[i, "CV"] <- CV
  
}

CV.pan$CV <- as.numeric(as.character(CV.pan$CV))

CV.pan %>% arrange(CV)

plot(CV.pan$CV)

gg.CV.pan <- CV.pan %>% mutate(color = ifelse(k == 1, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC", "Admixture.PAN.CV.pdf"), width = 4, height = 3.5)
gg.CV.pan
dev.off()

k <- 4
Q.k2.pan <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/", paste0("populations.26019snps.497ind.pan.recode.",2,".Q")))
Q.k3.pan <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/", paste0("populations.26019snps.497ind.pan.recode.",3,".Q")))
Q.k4.pan <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/", paste0("populations.26019snps.497ind.pan.recode.",4,".Q")))

Q.pan <- bind_rows(cbind(fam$V1, Q.k4.pan, K = 4),
                   cbind(fam$V1, Q.k3.pan, K = 3),
                   cbind(fam$V1, Q.k2.pan, K = 2))

head(Q.pan)

names(Q.pan) <- c("ID_GQ", paste0("Q", 1:4), "K")
#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str.pan <- Q.pan %>% pivot_longer(cols =  paste0("Q", 1:4), names_to = "Group", values_to = "Q") %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_col() +
  facet_grid(K ~ Region1, space = "free", scale = "free") +
  # scale_fill_manual(values = c("#00008bff","#436eeeff","#aa72dfff","#f4c9f4ff")) +
  scale_fill_manual(values = c("#00008bff","royalblue1","salmon","purple1")) +
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
gg.str.pan


# Figure: ordered by region

d.str.pan <- Q.pan %>% pivot_longer(cols =  paste0("Q", 1:4), names_to = "Group", values_to = "Q") %>%
  left_join(pop.data %>% select(ID_GQ, Region1, Location, Year, Month)) %>% 
  mutate(Region1 = factor(Region1, levels = c("REB","NHB","NWH","SWH","JAM","BEL","SEH","NEH","SHS","NHS","UNG","FRB","CSB","SLE")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>%
  mutate(Group.order = ifelse(is.na(Q), Group,
                              ifelse(K %in% 2, Group,
                                     ifelse(K %in% 3 & Group %in% "Q3", "Q1",
                                            ifelse(K %in% 3 & Group %in% "Q2", "Q3",
                                                   ifelse(K %in% 3 & Group %in% "Q1", "Q2",
                                                          ifelse(K %in% 4 & Group %in% "Q4", "Q1",
                                                                 ifelse(K %in% 4 & Group %in% "Q1", "Q4",
                                                                        Group))))))))

gg.str.pan <- d.str.pan %>% ggplot(aes(x = ID_GQ, y = Q, fill = Group.order)) + 
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(K ~ Region1, space = "free", scale = "free") +
  scale_fill_manual(values = c("#00008bff","royalblue1","salmon","purple1")) +
  # scale_fill_brewer(palette = "Set1") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        # axis.ticks.y = element_line(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC/", "Admixture_PAN_k2_to_k4.pdf"), width = 12, height = 6)
gg.str.pan
dev.off()

pop.data %>% subset(ID_GQ %in% Q.pan$ID_GQ) %>% left_join(d.str.pan) %>% subset(!duplicated(ID_GQ)) %>% group_by(Region1, Season) %>% summarise(N = n()) %>% print(n = 100)
pop.data %>% select(ID_GQ,Sex_qPCR,Region1,Location,Year,Month) %>% subset(ID_GQ %in% Q.pan$ID_GQ) %>% left_join(Q.pan) %>% subset(K %in% 3) %>%
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% arrange(Region1, Month) %>% View()
pop.data %>% subset(ID_GQ %in% Q.pan$ID_GQ) %>% left_join(Q.pan) %>% subset(Region1 %in% "SEH") %>% 
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))),
         Season = factor(Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))) %>% 
  group_by(Location, Season, Q2 > 0.7) %>% summarise(N = n()) %>% print(n = 100)

Q.pan.meta <- Q.pan %>% 
  left_join(pop.data %>% select(ID_GQ,Region1,Location,Year,Month,Day)) %>% 
  arrange(K, ID_GQ)
write.csv(Q.pan.meta, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC/Beluga_ddRAD_meta_admixture_PAN_k2-4.csv", row.names = F)


## Identify genetic clusters ----------------------------------------------

# K = 3
Q.pan 

hist(Q.pan$Q1[Q.pan$K %in% 3], breaks = 20)  # PAN cluster
table(Q.pan.meta$Region1[Q.pan.meta$Q1 > 0.5 & Q.pan.meta$K %in% 3])  # 7 CSB 15 FRB 2 JAM 15 NEH 20 NHB 24 NHS 21 NWH 25 SEH 82 SHS 2 SLE 11 SWH 70 UNG
table(Q.pan.meta$Region1[Q.pan.meta$Q1 > 0.5 & Q.pan.meta$K %in% 3 & Q.pan.meta$Month %in% c(7,8)])  # 5 CSB 8 FRB 15 NHB 4 NHS 13 NWH 22 SEH 22 SHS 1 SLE 11 SWH 50 UNG

hist(Q.pan$Q2[Q.pan$K %in% 3], breaks = 20)  # SEH cluster
table(Q.pan.meta$Region1[Q.pan.meta$Q2 > 0.5 & Q.pan.meta$K %in% 3])  # 20 SEH 1 SHS
table(Q.pan.meta$Region1[Q.pan.meta$Q2 > 0.5 & Q.pan.meta$K %in% 3 & Q.pan.meta$Month %in% c(7,8)])  # 17 SEH

hist(Q.pan$Q3[Q.pan$K %in% 3], breaks = 20)  # BEL-WHB cluster
table(Q.pan.meta$Region1[Q.pan.meta$Q3 > 0.5 & Q.pan.meta$K %in% 3])  # 22 BEL 3 JAM 1 NEH 16 NHB 2 NHS 40 NWH 5 SEH 12 HS 3 SWH 5 UNG
table(Q.pan.meta$Region1[Q.pan.meta$Q3 > 0.5 & Q.pan.meta$K %in% 3 & Q.pan.meta$Month %in% c(7,8)])  # 4 BEL 8 NHB 34 NWH 1 SEH 1 SHS 3 SWH 5 UNG

# Nither cluster is clearly defined

Q.pan %>% filter(K %in% 3) %>% nrow()  # 494 belugas
id.seh <- Q.pan %>% filter(Q2 > 0.5, K %in% 3) %>% pull(ID_GQ)
id.belwhb <- Q.pan %>% filter(Q3 > 0.5, K %in% 3) %>% pull(ID_GQ)
id.seh.belwhb <- c(id.seh, id.belwhb)  # 130
Q.pan %>% filter(K %in% 3, Q1 > 0.5) %>% nrow()  # 294 IDs, less than the 494-130 (364) belugas that aren't SEH or BEL-WHB

Q.pan.meta %>% filter(K %in% 3) %>% group_by(Region1) %>% summarise(Mean_Q1 = mean(Q1),
                                                                    Mean_Q2 = mean(Q2),
                                                                    Mean_Q3 = mean(Q3))
# Region1 Mean_Q1 Mean_Q2 Mean_Q3
# BEL       0.264  0.0940   0.642
# CSB       0.726  0.0950   0.179
# FRB       0.716  0.122    0.162
# JAM       0.315  0.0720   0.613
# NEH       0.574  0.176    0.250
# NHB       0.446  0.0668   0.487
# NHS       0.697  0.109    0.194
# NWH       0.348  0.0845   0.568
# SEH       0.367  0.387    0.247
# SHS       0.602  0.128    0.270
# SLE       0.629  0.163    0.208
# SWH       0.598  0.108    0.294
# UNG       0.616  0.119    0.265

Q.res.pan <- Q.pan.meta %>% filter(K %in% 3) %>% 
  mutate(Membership.pan = ifelse(Q2 > 0.5, "SEH",
                                 ifelse(Q3 > 0.5, "BEL-WHB",
                                        "PAN")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))))

write.csv(Q.res.pan, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02c_all_SNPs_HBSC/Beluga_ddRAD_meta_admixture_PAN_K3.csv", row.names = F)

## Fst, Ho and Fis analyses are not run for incipient clusters "identified" within th PAN cluster.
## Clusters within PAN are not clearly identified, and there is always a gradient for PAN.all, WHB.BEL, and SEH (albeit steep for SEH).
## TB: An important point is that this is what results look like at the whole genome level. So whatever this is, we aren't going to apply data at a higher
## resolution to resolve this. This is literally something super young or with gene flow between 'old' clusters.




# FST ---------------------------------------------------------------------
# FST within PAN clusters estimated using Regions instead of clusters, since they are of difficult definition: SEH vs WHB vs else

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
}

head(Q.res.pan)
Q.res.pan %>% group_by(Region1, Membership.pan) %>% summarise(N = n()) %>% print(n = 50)
Q.res.pan %>% group_by(Membership.pan) %>% summarise(N = n())
pop.admx.PAN.Q <- Q.res.pan %>% mutate(Region.pan = ifelse(Region1 %in% c("NHB","NWH","SWH","BEL"), "BEL-NWHB",
                                                           ifelse(Region1 %in% "SEH", "SEHB",
                                                                  "HS")))
pop.admx.PAN.Q %>% group_by(Region.pan, Region1) %>% summarise(N = n())

# SEH vs WHB vs else

## All IDs

gl.pan.3Q <- gl.pan.MAF05NA05
pop(gl.pan.3Q) <- data.frame(ID_GQ = indNames(gl.pan.3Q)) %>%
  left_join(pop.admx.PAN.Q) %>% pull(Region.pan)
table(pop(gl.pan.3Q), useNA = "ifany")

FST.3Q <- gl.fst.pop(gl.pan.3Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same

# ## 5 ID per cluster
# 
# id.fst.7Q <- c("S_20_01266","S_20_01428","S_20_01513","S_20_01530","S_20_04119",  # RES not present in ddRAD dataset
#                "S_20_01788","S_20_03480","S_20_01766","S_20_03664","S_20_01763",  # SLE
#                "S_20_00622","S_20_00710","S_20_00711","S_20_03448","S_20_03453",  # JAM
#                "S_20_00776","S_20_00661","S_20_00703","S_20_01418","S_20_03509",  # CSB
#                "S_20_01865","S_20_00621","S_20_01639","S_20_00990","S_20_03339",  # LGR
#                "S_20_00770","S_20_00830","S_20_00929","S_20_01063","S_20_02537",  # PANall (SWH,NWH,NHSmNHS,SHS)
#                "S_20_00678","S_20_00681","S_20_00809","S_20_01079","S_20_03515",  # PANwest (NWH,NWH,JAMmig,NWH,NWH)
#                "S_20_00817","S_20_01629","S_20_01708","S_20_03180","S_20_03658")  # SEH
# 
# length(id.fst.7Q)
# pop.admx.Q %>% subset(ID_GQ %in% id.fst.7Q) %>% nrow()
# pop.admx.Q %>% subset(ID_GQ %in% id.fst.7Q) %>% group_by(Region1) %>% summarise(N = n())
# pop.admx.Q %>% subset(ID_GQ %in% id.fst.7Q) %>% group_by(Membership) %>% summarise(N = n())
# 
# # d.fst.7Q <- pop.admx.Q %>% subset(ID_GQ %in% id.fst.7Q) %>% select(ID_GQ, Region1, Membership)
# # write.table(d.fst.7Q, "./02_Results/01_ddRAD_Bringloe/03_FST/Fst_pops_ddRAD_lcWGS_230518.txt")
# 
# gl.red.7Q <- gl.keep.ind(gl.final.7Q, id.fst.7Q, recalc = T, verbose = 5)  # keep only 35 individual of choice
# table(pop(gl.red.7Q), useNA = "ifany")
# FST.red.7Q <- gl.fst.pop(gl.red.7Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same
# 
# save(list = c("FST.3Q"),
#      file = file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_pan.Rdata"))
load(file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_pan.Rdata"))

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


## Figures: SLE vs CSB vs JAM-BEL vs LGR vs PANall vs PANwest vs SEH ------

### ddRAD -----------------------------------------------------------------

FST.table.3Q <- FST.3Q %>% table.fst()
# Population1 Population2 Lower bound CI limit Upper bound CI limit p-value         Fst
#    BEL-NWHB          HS         0.0009781981          0.001119884       0 0.001045126
#    BEL-NWHB        SEHB         0.0024928526          0.002786749       0 0.002639936
#          HS        SEHB         0.0018724994          0.002103334       0 0.001984748


FST.3Q %>% table.fst() %>%
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst),
            Max.Pvalue = round(max(`p-value`),))
#        Mean           sd         Min         Max Max.Pvalue
# 0.001889937 0.0008016209 0.001045126 0.002639936          0

FST.ddRAD.3Q <- FST.3Q %>% heat.fst() %>%
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)

gFST.ddRAD.3Q <-  FST.ddRAD.3Q %>% 
  mutate(Population1 = factor(Population1, levels = c("BEL-NWHB","SEHB","HS")),
         Population2 = factor(Population2, levels = c("BEL-NWHB","SEHB","HS"))) %>%
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) +
  geom_tile(colour=  "white") +
  geom_text(aes(label = round(Fst, digits = 3)), color = "black", size = 9) +  # fontface = "bold"
  scale_fill_viridis_c(name = "Fst", na.value = "white")+
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.key.size = unit(1.25, "cm"),
        legend.title = element_text(size = 30),  # , face = "bold"
        legend.text = element_text(size = 22)) +
  theme(strip.text = element_text(angle = 0),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 27.5, color = "black"))
# pdf(file = "02_Results/01_ddRAD_Bringloe/03_FST/03c_all_SNPs_HBSC/FST.3Q.ddRAD.HBSC.230725.pdf", width = 9, height = 7.5)
# png(file = "02_Results/01_ddRAD_Bringloe/03_FST/03c_all_SNPs_HBSC/FST.3Q.ddRAD.HBSC.230725.png", width = 9, height = 7.5, units = "in", res = 150, pointsize = 12)
gFST.ddRAD.3Q
dev.off()



