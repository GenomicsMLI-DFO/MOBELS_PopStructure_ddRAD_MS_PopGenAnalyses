# Info --------------------------------------------------------------------
#
# Overview: Population genetics analyses for eastern Arctic beluga Delphinapterus leucas using ddRAD loci
# 
# Authors: Luca Montana and Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2022-11-02
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

## What is the sample size
pop.data %>% pull(Cat_sample) %>% table()
pop.data %>% pull(ID) %>% table()
pop.data %>% pull(Region) %>% table(useNA = "ifany")  # all regions
pop.data %>% pull(Region1) %>% table(useNA = "ifany")  # JAM and LON merged into JAM
pop.data %>% pull(Common_name) %>% table()

pop.data %>% filter(Cat_sample %in% "Sample") %>% 
  group_by(Region2) %>%  # NWH, SWH, NHB and REB merged into WHB
  summarise(N = n()) #%>% write_csv("clipboard")

pop.data %>% filter(Cat_sample %in% "Sample") %>% 
  group_by(Month) %>% 
  summarise(N = n()) #%>% write_csv("clipboard")


## ddRAD ------------------------------------------------------------------

if(!file.exists(file.path("./00_Data", "03_ddRAD_Bringloe"))){
  dir.create(file.path("./00_Data", "03_ddRAD_Bringloe"))
  print(file.path("./00_Data", "03_ddRAD_Bringloe"))
}

load("./00_Data/01_Filtering.ref/01_Bringloe/01i_UniqueFinal/populations.60102snps.638ind.H06.DP.single.final.recode.vcf.adegenet.Rdata")
# load("../MOBELS_PopStructure_ddRAD_MS_SNPdiscovery_noproj/00_Data/06b_Filtering.ref/1_Beluga/06g_UniqueFinal/populations.59156snps.650ind.H06.DP.single.final.recode.vcf.adegenet.Rdata")

# gl.data <- gl.final
# gi.data <- gi.final  # use this if you want to exclude SLE

pop(gl.final) <- data.frame(ID_GQ = indNames(gl.final)) %>% 
  left_join(pop.data) %>% pull(Region1)
table(pop(gl.final), useNA = "ifany")
# BEL CSB FRB JAM NEH NHB NHS NWH SEH SHS SLE SWH UNG 
#  43  27  16  24  32  39  27  68 124 112  23  14  89

pop(gi.final) <- data.frame(ID_GQ = indNames(gi.final)) %>% 
  left_join(pop.data) %>% pull(Region1)
table(pop(gi.final), useNA = "ifany")

# Subset Minor Allele Frequency
maf <- c(0.01,0.05,0.1)

## All seasons
LOC.MAF.NA <- lapply(maf, function(i){
  loc.maf.na <- filter.MAF.NA(gi.final, MAF.trs = i, NA.trs = 0.05)
  loc.maf.na
})

LOC.MAF01.NA05 <- LOC.MAF.NA[[1]]  # 45571 loci
LOC.MAF05.NA05 <- LOC.MAF.NA[[2]]  # 26019 loci
LOC.MAF10.NA05 <- LOC.MAF.NA[[3]]  # 15675 loci

gl.final.MAF01NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF01.NA05]
gi.final.MAF01NA05 <- gi.final[loc = LOC.MAF01.NA05]
hist(adegenet::minorAllele(gi.final.MAF01NA05))
locMAF01NA05 <- data.frame("ID" = locNames(gl.final.MAF01NA05))  # 45571 loci

gl.final.MAF05NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF05.NA05]
gi.final.MAF05NA05 <- gi.final[loc = LOC.MAF05.NA05]
hist(adegenet::minorAllele(gi.final.MAF05NA05))
locMAF05NA05 <- data.frame("ID" = locNames(gl.final.MAF05NA05))  # 26019 loci

gl.final.MAF10NA05 <- gl.final[, locNames(gl.final) %in% LOC.MAF10.NA05]
gi.final.MAF10NA05 <- gi.final[loc = LOC.MAF10.NA05]
hist(adegenet::minorAllele(gi.final.MAF10NA05))
locMAF10NA05 <- data.frame("ID" = locNames(gl.final.MAF10NA05))  # 15675 loci

# Save loci list to filter vcf file later

write.csv(locMAF01NA05, file.path("./00_Data/03_ddRAD_Bringloe", "Loc.MAF01NA05.csv"), row.names = F, quote = F)
write.csv(locMAF05NA05, file.path("./00_Data/03_ddRAD_Bringloe", "Loc.MAF05NA05.csv"), row.names = F, quote = F)
write.csv(locMAF10NA05, file.path("./00_Data/03_ddRAD_Bringloe", "Loc.MAF10NA05.csv"), row.names = F, quote = F)




# Stats on sampling size --------------------------------------------------

pop.data %>% filter(ID_GQ %in% indNames(gi.final)) %>% group_by(Region, Month) %>% summarise(N = n()) %>% print(n = 80)
# CSB   Jun (1) Jul (21) Aug (4) NA (1)
# FRB   Jun (8) Jul (7) Aug (1)
# JAM   Jul (3) Aug (15) Sep (5)
# LON   Oct (1)
# NEH   May (2) Sep (1) Oct (28) NA (1)
# NHB   Jul (7) Aug (17) Sep (13) NA (2)
# NHS   May (1) Jun (3) Jul (5) Oct (12) Nov (6)
# NWH   Jul (10) Aug (42) Sep (13) NA (3)
# REB   Aug (5) Sep (7) NA (1)
# SAN   May (21) Jun (8) Jul (5) Sep (5) Dec (2) NA (2)
# SEH   Jun (12) Jul (46) Aug (54) Sep (2) Oct (1) NA (9)
# SHS   Jun (43) Jul (23) Oct (40) Nov (4) NA (2)
# SLE   May (3)  Jun (1)  Jul (11) Aug (3) Sep (2) Oct (2) Nov (1)
# SWH   Jul (11) Aug (3)
# UNG   Jun (21) Jul (61) Aug (6) Sep (1)

## All seasons

stat.gen <- data.frame(ID_GQ = indNames(gl.final)) %>% left_join(pop.data)

stat.gen %>% group_by(Year) %>% summarise(N = n()) %>% print(n = 50)  # 1989-2011, 2015-2019
stat.gen %>% group_by(Region1) %>% summarise(N = n(), Min_year = min(Year), Max_year = max(Year))
# Region1     N Min_year Max_year
# BEL        43     2002     2005
# CSB        27     2002     2007
# FRB        16     1991     1992
# JAM        24     2002     2009
# NEH        32     1998     2018
# NHB        39     1993     2006
# NHS        27     1989     2000
# NWH        68     1992     2015
# SEH       124     1990     2018
# SHS       112     1994     2015
# SLE        23     2000     2019
# SWH        14     2002     2005
# UNG        89     1994     2018
stat.gen %>% filter(Year %in% NA) %>% View()  # all years are known

stat.gen %>% group_by(Region1, Month) %>% summarise(N = n()) %>% print(n = 100)

stat.gen %>% group_by(Region) %>% summarise(N = n())
stat.gen %>% group_by(Region1) %>% summarise(N = n())
stat.gen %>% group_by(Region2) %>% summarise(N = n())
# Region2     N
# BEL        43
# CSB        27
# FRB        16
# JAM        24
# NEH        32
# NHS        27
# SEH       124
# SHS       112
# SLE        23
# UNG        89
# WHB       121

stat.gen %>% group_by(Region2) %>% summarise(N = n_distinct(Haplo))

head(stat.gen)

summary(stat.gen)




# Filter ddRAD dataset ----------------------------------------------------

# Filter by specimens
ind.names <- sort(indNames(gi.final))
ind.names  # N = 638


# MAF and NA

## MAF01NA05 (locMAF01NA05 - about 45k loci)

cmd1a <- paste("--vcf", file.path("./00_Data/01_Filtering.ref/01_Bringloe/01i_UniqueFinal/populations.60102snps.638ind.single.final.recode.vcf"),
               "--recode",
               "--snps", file.path("./00_Data/03_ddRAD_Bringloe/", "Loc.MAF01NA05.csv"),
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", paste0("populations.",nrow(locMAF01NA05),"snps.",length(ind.names),"ind.final"))
)
cmd1a

A1a <- system2("vcftools", cmd1a, stdout=T, stderr=T)
tail(A1a)


## MAF05NA05 (locMAF05NA05 - about 26k loci)

cmd1b <- paste("--vcf", file.path("./00_Data/01_Filtering.ref/01_Bringloe/01i_UniqueFinal/populations.60102snps.638ind.single.final.recode.vcf"),
               "--recode",
               "--snps", file.path("./00_Data/03_ddRAD_Bringloe", "Loc.MAF05NA05.csv"),
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", paste0("populations.",nrow(locMAF05NA05),"snps.",length(ind.names),"ind.final"))
)
cmd1b

A1b <- system2("vcftools", cmd1b, stdout=T, stderr=T)
tail(A1b)


## MAF10NA05 (locMAF10NA05 - about 15k loci)

cmd1c <- paste("--vcf", file.path("./00_Data/01_Filtering.ref/01_Bringloe/01i_UniqueFinal/populations.60102snps.638ind.single.final.recode.vcf"),
              "--recode",
               "--snps", file.path("./00_Data/03_ddRAD_Bringloe", "Loc.MAF10NA05.csv"),
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", paste0("populations.",nrow(locMAF10NA05),"snps.",length(ind.names),"ind.final"))
)
cmd1c

A1c <- system2("vcftools", cmd1c, stdout=T, stderr=T)
tail(A1c)


## Save as plink (admixture) tped (for pcadapt) ---------------------------

# MAF05NA05

cmd2b <- paste("--vcf", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.vcf"), 
               #"--recode",
               "--plink-tped",
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.638ind.final.recode")
)
cmd2b

A2b <- system2("vcftools", cmd2b, stdout=T, stderr=T)
A2b


## Make bed files ---------------------------------------------------------

# MAF05NA05

cmd3b <- paste("--tfam", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.tfam"), 
               "--tped", file.path("./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.tped"), 
               "--make-bed", 
               "--out", file.path("./00_Data/03_ddRAD_Bringloe", "populations.26019snps.638ind.final.recode")
               
)
cmd3b

A3b <- system2(file.path(plink_path, "plink"), cmd3b, stdout=T, stderr=T)
A3b

pop.gen.meta <- pop.data %>% subset(ID_GQ %in% ind.names) %>% 
  droplevels() %>%  # to remove unused Pop levels
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% 
  select(ID_GQ,ID,ID_rcpt,Cat_sample,No_plate,No_well,Barcode,Genus,Species,Age,Sex_visual,Sex_qPCR,Region,Region1,Region2,Location,Lat,Lon,Year,Month,Day,Community,
         Haplo,Season)

write.csv(pop.gen.meta, file = "./00_Data/03_ddRAD_Bringloe/Beluga_ddRAD_analyses.csv", row.names = F)




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

l.gl <- lapply(ls(pattern = "gl.final.MAF05NA05"), function(x) get(x))  # put genlight object in a list to perform PCA using a lapply
names(l.gl) <- ls(pattern = "gl.final.MAF05NA05")  # name genlight objects within list

# Prepare GI objects: useful for estimation of sample heterozygosity

## All regions - all seasons

l.gi <- lapply(ls(pattern = "gi.final.MAF05NA05"), function(x) get(x))  # put genlight object in a list to perform PCA using a lapply
names(l.gi) <- ls(pattern = "gi.final.MAF05NA05")  # name genlight objects within list


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

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01a_all_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01a_all_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/01_PCA", "01a_all_SNPs"))
}


# save(list = c("l.pca"), file = "./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_all.Rdata")
load("./02_Results/01_ddRAD_Bringloe/01_PCA/PCA_all.Rdata")

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/", ".gitignore")) ){
  cat("*.Rdata", "*.RData", "!.gitignore", sep = "\n",
      file = file.path("./02_Results/01_ddRAD_Bringloe/01_PCA/", ".gitignore")) 
}

# Eig var

MAF05NA05.var.638 <- pca_var(l.pca[["gl.final.MAF05NA05"]], nInd(gl.final)-1) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.3,0.6,0.9,1.2)) +
  # geom_hline(yintercept = 0, col = "grey20") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
MAF05NA05.var.638

MAF05NA05.var.30 <- pca_var(l.pca[["gl.final.MAF05NA05"]], 30) %>% 
  ggplot(aes(x = axis, y = p.eig * 100, fill = axis)) +
  plot_pca_eig +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.3,0.6,0.9,1.2)) +
  # geom_hline(yintercept = 0, col = "grey20") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
MAF05NA05.var.30

# With inset
MAF05NA05.var <- MAF05NA05.var.30 + annotation_custom(ggplotGrob(MAF05NA05.var.638), xmin = 13, xmax = 30,
                                     ymin = 0.6, ymax = 1.18)
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/Var.PCA.global.inset.MAF05NA05.230717.pdf", width = 13, height = 10)
MAF05NA05.var
dev.off()



## Figures (MAF05NA05) ----------------------------------------------------

region.labels <- c("Belcher Islands","Cumberland Sound","Frobisher Bay","South-East Hudson Bay","James Bay","North-East Hudson Bay","North Hudson Strait",
                   "South Hudson Strait","Saint Lawrence estuary","Ungava Bay","Western Hudson Bay")
names(region.labels) <- c("BEL","CSB","FRB","SEH","JAM","NEH","NHS","SHS","SLE","UNG","WHB")

# All regions and seasons

all.pca <- l.pca[["gl.final.MAF05NA05"]] %>% QuickPop::pca_scoretable(naxe = 8) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>%
  droplevels() %>%  # to remove unused Pop levels
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))))

## All-in-one

table(all.pca$Region2)
all.pca$Region2 <- factor(all.pca$Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","SEH","BEL","JAM","WHB"))
table(all.pca$Season, useNA = "ifany")
all.pca$Season <- factor(all.pca$Season, levels = c("Spring","Summer","Fall","Winter","Unknown"))

gPCA.MAF05NA05.1v2 <- ggplot(data = all.pca, aes(x = score.PC1, y = score.PC2, col = Region2, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1) +  # added stroke for presentation
  # scale_color_manual(name = "Harvest region", values = c(alpha("#79c9e4ff",0.75),alpha("#49c687ff",0.75),alpha("#f4c9f4ff",0.75),alpha("#aa72dfff",0.75),
  #                                                        alpha("#B75BA7",0.75),alpha("#C64D83",0.75),alpha("#EBA783",0.75),alpha("#f6835f",0.75),
  #                                                        alpha("#00008bff",0.75),alpha("#e5b65fff",0.75),alpha("#436eeeff",0.75)),
  #                    labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
  #                               "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_color_manual(name = "Harvest region", values = c(alpha("deepskyblue",0.65),alpha("springgreen4",0.65),alpha("plum2",0.65),alpha("purple3",0.65),
                                                         alpha("orchid3",0.65),alpha("palevioletred2",0.65),alpha("salmon",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("orange",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none",
         size = "none",
         shape = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.MAF05NA05.1v2.20230717.pdf", width = 13, height = 10)
gPCA.MAF05NA05.1v2
dev.off()

gPCA.MAF05NA05.3v4 <- ggplot(data = all.pca, aes(x = score.PC3, y = score.PC4, col = Region2, shape = Season, size = Season)) +
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
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none",
         size = "none",
         shape = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.MAF05NA05.3v4.20230717.pdf", width = 13, height = 10)
gPCA.MAF05NA05.3v4
dev.off()
# du.pca %>% filter(Pop %in% "JAM") %>% group_by(score.PC4<(-10), Location) %>% summarise(N = n()) %>% print(n = 100)
# # 4 JAM specimens with PC4 < -10. From Pointe de Repentigny (but other from Pointe de Repentigny in main cluster)
# du.pca %>% filter(Pop %in% "JAM") %>% group_by(score.PC4<(-10), Location, Year) %>% summarise(N = n()) %>% print(n = 100)
# # All in 2007 at Pointe de Repentigny (midway between Eastmain and Wemindji), but another speicment from main cluster harvested there in 2007


gPCA.MAF05NA05.5v6 <- ggplot(data = all.pca, aes(x = score.PC5, y = score.PC6, col = Region2, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(stroke = 1) +  # added stroke for presentation
  scale_y_continuous(breaks = c(-15,-10,-5,0,5,10)) +
  scale_x_continuous(breaks = c(-15,-10,-5,0,5,10,15)) +
  scale_color_manual(name = "Harvest region", values = c(alpha("deepskyblue",0.65),alpha("springgreen4",0.65),alpha("plum2",0.65),alpha("purple3",0.65),
                                                         alpha("orchid3",0.65),alpha("palevioletred2",0.65),alpha("salmon",0.65),alpha("red3",0.65),
                                                         alpha("chocolate3",0.65),alpha("orange",0.65),alpha("royalblue1",0.65)),
                     labels = c("Saint Lawrence estuary","Cumberland Sound","Frobisher Bay","Ungava Bay","South Hudson Strait","North Hudson Strait",
                                "North-East Hudson Bay","South-East Hudson Bay","Belcher Islands","James Bay","Western Hudson Bay")) +
  scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # squares for spring-fall, 16 for summer, triangles for NA
  labs(x = paste0("PC5 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[5] %>% round(3) * 100, "%)"),
       y = paste0("PC6 (", QuickPop::pca_var(l.pca[["gl.final.MAF05NA05"]])$p.eig[6] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none",
         size = "none",
         shape = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.MAF05NA05.5v6.20230717.pdf", width = 13, height = 10)
gPCA.MAF05NA05.5v6
dev.off()
# du.pca %>% filter(Pop %in% "EHB") %>% group_by(score.PC5>4, Location) %>% summarise(N = n()) %>% print(n = 100)
# # 24 EHB specimens with PC5 > 4. No clear spatial pattern
# du.pca %>% filter(Pop %in% "EHB") %>% group_by(score.PC5>4, Location, Year) %>% summarise(N = n()) %>% print(n = 100)
# # Nor spatio-temporal pattern to distinguish them from other EHB specimens


## Figures (MAF05NA05) - NA -----------------------------------------------

l.na.info <- lapply(l.gl, function(w){
  na.info <- data.frame(ID_GQ = indNames(w),
                        NNA = count.ind.na.gl(w))
  na.info
})


gPCA.NA.MAF05NA05.1v2 <- all.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, color = NNA)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[[1]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[[1]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_colourbar(title = "% NA", order = 1)) +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(legend.position = c(0.825,0.800),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.1v2

gPCA.NA.MAF05NA05.3v4 <- all.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, color = NNA)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[[1]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[[1]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.3v4

gPCA.NA.MAF05NA05.5v6 <- all.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6, color = NNA)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5) +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC5 (", QuickPop::pca_var(l.pca[[1]])$p.eig[5] %>% round(3) * 100, "%)"),
       y = paste0("PC6 (", QuickPop::pca_var(l.pca[[1]])$p.eig[6] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.NA.MAF05NA05.5v6

gPCA.NA.1v2.3v4.5v6 <- ggpubr::ggarrange(gPCA.NA.MAF05NA05.1v2,
                                         gPCA.NA.MAF05NA05.3v4,
                                         gPCA.NA.MAF05NA05.5v6,
                                         nrow = 3, ncol = 1,
                                         common.legend = T, align = "hv", legend = "right")
gPCA.NA.1v2.3v4.5v6




## Figures (MAF05NA05) - Heterozygosity -----------------------------------

l.ho <- lapply(l.gl, function(h){
  ho <- gl.report.heterozygosity(h, method='ind')  # Verify what gl.report.heterozigosity does
  ho
})


gPCA.Ho.MAF05NA05.1v2 <- all.pca %>% 
  left_join(l.ho[[1]], by = c("ID" = "ind.name")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, color = Ho)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC1 (", QuickPop::pca_var(l.pca[[1]])$p.eig[1] %>% round(3) * 100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(l.pca[[1]])$p.eig[2] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = guide_colourbar(title = "Observed\nheterozygosity", order = 1)) +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(legend.position = c(0.825,0.800),
        legend.spacing.y = unit(0.25, "cm"),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "horizontal") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.1v2

gPCA.Ho.MAF05NA05.3v4 <- all.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, color = NNA)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC3 (", QuickPop::pca_var(l.pca[[1]])$p.eig[3] %>% round(3) * 100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(l.pca[[1]])$p.eig[4] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.3v4


gPCA.Ho.MAF05NA05.5v6 <- all.pca %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6, color = NNA)) +
  # ggplot(aes(x = score.PC1, y = score.PC2, color = NNA, shape = Season, size = Season)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  # geom_point(stroke = 1) +
  geom_point(size = 7, stroke = 1) +
  viridis::scale_color_viridis(alpha = 0.5, option = "plasma") +
  # scale_size_manual(values = c(3,7.5,3,3,3)) +  # 3 for spring-fall and NA, 9 for summer
  scale_shape_manual(values = c(6,19,2,0,5)) +  # reversed triangle Spring, dots Summer, tringle Fall, squadre Winter, losange Unknown
  # facet_wrap(~ Region2, labeller = labeller(Region2 = region.labels)) +
  labs(x = paste0("PC5 (", QuickPop::pca_var(l.pca[[1]])$p.eig[5] %>% round(3) * 100, "%)"),
       y = paste0("PC6 (", QuickPop::pca_var(l.pca[[1]])$p.eig[6] %>% round(3) * 100, "%)")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(colour = "none") +
  # guides(colour = guide_colourbar(title = "% NA", order = 1),
  #        shape = guide_legend(title = "Harvest season", order = 2),
  #        size = "none") +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gPCA.Ho.MAF05NA05.5v6

gPCA.NNA.Ho.1v2.3v4.5v6 <- ggpubr::ggarrange(gPCA.NA.MAF05NA05.1v2,gPCA.Ho.MAF05NA05.1v2,
                                             gPCA.NA.MAF05NA05.3v4,gPCA.Ho.MAF05NA05.3v4,
                                             gPCA.NA.MAF05NA05.5v6,gPCA.Ho.MAF05NA05.5v6,
                                             nrow = 3, ncol = 2,
                                             common.legend = F)
# pdf(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.MAF05NA05.NNA.HO.1v2.3v4.5v6.20230717.pdf", width = 15, height = 20)
# png(file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.MAF05NA05.NNA.HO.1v2.3v4.5v6.20230717.png", width = 15, height = 20, units = "in", res = 150, pointsize = 12)
gPCA.NNA.Ho.1v2.3v4.5v6
dev.off()



## Save results PCA -------------------------------------------------------

PCA.all.MAF05NA05 <- all.pca %>%
  left_join(l.ho[[1]], by = c("ID" = "ind.name")) %>% 
  left_join(l.na.info[[1]], by = c("ID" = "ID_GQ")) %>%
  mutate(Sex_qPCR = factor(Sex_qPCR, levels = c("F","M","NA")))

write.csv(PCA.all.MAF05NA05, file = "02_Results/01_ddRAD_Bringloe/01_PCA/01a_all_SNPs/PCA.all.MAF05NA05_n638.csv", row.names = F)




# Admixture ---------------------------------------------------------------

## All Pops ---------------------------------------------------------------

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)

# Admixture analysis

# set.seed(111)

for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("Bringloe.k",k, ".log"),
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
  temp <- readLines(file.path("./02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("Bringloe.k",k, ".log")))
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
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", "Admixture.all.CV.pdf"), width = 4, height = 3.5)
# png(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", "Admixture.all.CV.png"), width = 4, height = 3.5, units = "in", res = 150, pointsize = 12)
gg.CV
dev.off()

k <- 6
Q.k2.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",2,".Q")))
Q.k3.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",3,".Q")))
Q.k4.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",4,".Q")))
Q.k5.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",5,".Q")))
Q.k6.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",6,".Q")))
Q.k7.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",7,".Q")))
Q.k8.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",8,".Q")))
Q.k9.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",9,".Q")))
Q.k10.all <-  read.table(file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", paste0("populations.26019snps.638ind.final.recode.",10,".Q")))


Q.all <- bind_rows(cbind(fam$V1, Q.k10.all, K = 10),
                   cbind(fam$V1, Q.k9.all, K = 9),
                   cbind(fam$V1, Q.k8.all, K = 8),
                   cbind(fam$V1, Q.k7.all, K = 7),
                   cbind(fam$V1, Q.k6.all, K = 6),
                   cbind(fam$V1, Q.k5.all, K = 5),
                   cbind(fam$V1, Q.k4.all, K = 4),
                   cbind(fam$V1, Q.k3.all, K = 3),
                   cbind(fam$V1, Q.k2.all, K = 2))

head(Q.all)
names(Q.all) <- c("ID_GQ", paste0("Q", 1:10), "K")

# Figure (check this out: https://www.royfrancis.com/pophelper/articles/index.html)

gg.str.all <- Q.all %>% pivot_longer(cols =  paste0("Q", 1:10), names_to = "Group", values_to = "Q") %>% 
  left_join(pop.data) %>% 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  geom_col() +
  facet_grid(K ~ Region2 , space = "free", scale = "free") +
  scale_fill_brewer(palette = "Set1") +
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


# Identify clusters

Q.all <- Q.all %>%
  left_join(pop.data %>% select(ID_GQ,ID_rcpt,Region1,Location,Year,Month,Day,Community,Haplo)) %>%
  arrange(K, ID_GQ)
write.csv(Q.all, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/Beluga_ddRAD_meta_admixture_ALL_k2-10.csv", row.names = F)

## K = 5
hist(Q.all$Q1[Q.all$K %in% 5], breaks = 20)  # SLE very differentiated (Q > 0.9)

hist(Q.all$Q3[Q.all$K %in% 5], breaks = 20)  # CSB
table(Q.all$Region1[Q.all$Q3 > 0.7 & Q.all$K %in% 5])
Q.all %>% subset(K %in% 5 & Region1 %in% "CSB") %>% View()
Q.all %>% subset(K %in% 5) %>% View()  # 20 CSB, 1 WHB, 1 FRB (consistent with ResDoc), possible hybrids from CSB and WHB

hist(Q.all$Q4[Q.all$K %in% 5], breaks = 20)  # JAM-BEL very differentiated
Q.all %>% subset(K %in% 5) %>% View()  # Q4 > 0.75 for clear JAM-BEL identification
Q.all %>% subset(K %in% 5 & Region1 %in% "JAM") %>% View()  # All JAM Q4 > 0.9 and LON Q4 = 0.83
Q.all %>% subset(K %in% 5 & Region1 %in% "BEL") %>% View()  # 4 SAN Q4 > 0.9 (+ 1 at 0.89). Overall gradient from 0.75 to 0.95. Difficult to distinguish from JAM

hist(Q.all$Q5[Q.all$K %in% 5], breaks = 20)  # LGR very differentiated (and 2 hybrids Q5 = 0.36-37)
table(Q.all$Region1[Q.all$Q5 > 0.70 & Q.all$K %in% 5])  # 34 EHB, 11 NEH, 1 NHS, 11 SHS (0.74: 34 EHB, 10 NEH, 1 NHS, 7 SHS, no WHB; 0.6: 34 EHB, 11 NEH, 1 NHS, 12 SHS, 1 WHB)
Q.all %>% subset(K %in% 5) %>% View()  # A few individuals from NEH, SHS, and WHB betwee 0.6 and 0.74 ( I consider those above 0.74 as 'real' LGR individuals harvested during migrations)
Q.all %>% subset(K %in% 5 & Region1 %in% "EHB") %>% View()  #Q5 > 0.74 for clear LGR identification.

hist(Q.all$Q2[Q.all$K %in% 5], breaks = 20)  # LGR very differentiated (and 2 hybrids Q5 = 0.36-37)

Q.all.K5 <- Q.all %>% subset(K %in% 5) %>% 
  mutate(Membership = ifelse(Q1 > 0.7, "SLE",
                             ifelse(Q3 > 0.7, "CSB",
                                    ifelse(Q4 > 0.7, "JAM",
                                           ifelse(Q5 > 0.7, "LGR",
                                                  "PAN")))))


# Figure: ordered by genetic cluster and region

d.str.all <- Q.all %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>%
  left_join(Q.all.K5 %>% select(ID_GQ, Membership)) %>% 
  mutate(Region1 = factor(Region1, levels = c("NHB","NWH","SWH","JAM","BEL","SEH","NEH","SHS","NHS","UNG","FRB","CSB","SLE")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))),
         Membership = factor(Membership, levels = c("PAN","LGR","CSB","JAM","SLE"))) %>% 
  mutate(Group.order = ifelse(is.na(Q), Group,
                         ifelse(K %in% 2, Group,
                                ifelse(K %in% 3 & Group %in% "Q1", "Q3",
                                ifelse(K %in% 3 & Group %in% "Q3", "Q1",
                                       ifelse(K %in% 4 & Group %in% "Q1", "Q3",
                                       ifelse(K %in% 4 & Group %in% "Q3", "Q2",
                                       ifelse(K %in% 4 & Group %in% "Q2", "Q1",
                                              ifelse(K %in% 5 & Group %in% "Q1", "Q2",
                                              ifelse(K %in% 5 & Group %in% "Q3", "Q5",
                                              ifelse(K %in% 5 & Group %in% "Q4", "Q3",
                                              ifelse(K %in% 5 & Group %in% "Q5", "Q4",
                                              ifelse(K %in% 5 & Group %in% "Q2", "Q1",
                                                     ifelse(K %in% 6 & Group %in% "Q6", "Q1",
                                                     ifelse(K %in% 6 & Group %in% "Q4", "Q3",
                                                     ifelse(K %in% 6 & Group %in% "Q3", "Q4",
                                                     ifelse(K %in% 6 & Group %in% "Q1", "Q6",
                                                            Group)))))))))))))))))

gg.str.all <- d.str.all %>% filter(K %in% c(2:6)) %>% ggplot(aes(x = ID_GQ, y = Q, fill = Group.order)) + 
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(K ~ Membership + Region1 , space = "free", scale = "free") +
  # scale_fill_manual(values=c("#436eeeff", "#79c9e4ff", "#e5b65fff", "#f6835f", "#49c687ff", "#00008bff")) +
  scale_fill_manual(values=c("royalblue1", "deepskyblue", "orange", "red3", "springgreen4", "darkblue")) +
  # scale_fill_brewer(palette = "Set1") +
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
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"))
# pdf(file = file.path(here::here(), "02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/", "Admixture_k2_to_k6.pdf"), width = 12, height = 6)
gg.str.all
dev.off()


## Identify genetic clusters ----------------------------------------------

# K = 5
Q.all 

## SLE
hist(Q.all$Q1[Q.all$K %in% 5], breaks = 20)  # SLE very differentiated (Q > 0.9)

## CSB
hist(Q.all$Q3[Q.all$K %in% 5], breaks = 20)  # CSB
table(Q.all$Region1[Q.all$Q3 > 0.7 & Q.all$K %in% 5])    # 20 CSB, 1 WHB, 1 FRB (consistent with ResDoc), possible hybrids from CSB and WHB
Q.all %>% subset(K %in% 5 & Region1 %in% "CSB") %>% View()
Q.all %>% subset(K %in% 5) %>% View()

## JAM
hist(Q.all$Q4[Q.all$K %in% 5], breaks = 20)  # JAM very differentiated
Q.all %>% subset(K %in% 5) %>% View()  # Q4 > 0.75 for clear JAM identification
Q.all %>% subset(K %in% 5 & Region1 %in% "JAM") %>% View()  # All JAM Q4 > 0.9 and LON Q4 = 0.83
Q.all %>% subset(K %in% 5 & Region1 %in% "BEL") %>% View()  # 4 SAN Q4 > 0.9 (+ 1 at 0.89). Overall gradient from 0.75 to 0.95. Difficult to distinguish from JAM

## LGR
hist(Q.all$Q5[Q.all$K %in% 5], breaks = 20)  # LGR very differentiated (and 2 hybrids Q5 = 0.36-37)
table(Q.all$Region1[Q.all$Q5 > 0.70 & Q.all$K %in% 5])  # 34 SEH, 11 NEH, 1 NHS, 11 SHS
Q.all %>% subset(K %in% 5) %>% View()  # A few individuals from NEH, SHS, and WHB between 0.6 and 0.74 ( I consider those above 0.74 as 'real' LGR individuals harvested during migrations)
Q.all %>% subset(K %in% 5 & Region1 %in% "SEH") %>% View()  #Q5 > 0.74 for clear LGR identification; 0.7 is a good compromise

## PAN
hist(Q.all$Q2[Q.all$K %in% 5], breaks = 20)  # PAN

Q.res <- Q.all %>% subset(K %in% 5) %>% 
  mutate(Membership.all = ifelse(Q1 > 0.5, "SLE",
                                 ifelse(Q3 > 0.5, "CSB",
                                        ifelse(Q4 > 0.5, "JAM",
                                               ifelse(Q5 > 0.5, "LGR",
                                                      "PAN"))))) %>% 
  select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Community, Haplo, K, Q1, Q2, Q3, Q4, Q5, Membership.all)

write.csv(Q.res, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/02a_all_SNPs/Beluga_ddRAD_meta_admixture_ALL_K5.csv", row.names = F)




# FST ---------------------------------------------------------------------

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "03_FST"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03a_all_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03a_all_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe/03_FST", "03a_all_SNPs"))
}


head(Q.res)
Q.res %>% group_by(Membership.all) %>% summarise(N = n())

# SLE vs CSB vs JAM vs LGR vs PAN

## All IDs

gl.final.5Q <- gl.final.MAF05NA05
pop(gl.final.5Q) <- data.frame(ID_GQ = indNames(gl.final.5Q)) %>% 
  left_join(Q.res) %>% pull(Membership.all)
table(pop(gl.final.5Q), useNA = "ifany")

FST.5Q <- gl.fst.pop(gl.final.5Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same

## 5 ID per cluster: all in common between lcWGS and ddRAD datasets but for LGR cluster (only two in common)

id.fst.5Q <- c(#"S_20_01266","S_20_01428","S_20_01513","S_20_01530","S_20_04119",  # RES not present in ddRAD dataset
  "S_20_01788","S_20_03480","S_20_01766","S_20_03664","S_20_01763",  # SLE
  "S_20_00622","S_20_00710","S_20_00711","S_20_03448","S_20_03453",  # JAM
  "S_20_00776","S_20_00661","S_20_00703","S_20_01418","S_20_03509",  # CSB
  "S_20_01865","S_20_00621","S_20_01639","S_20_00990","S_20_03339",  # LGR
  "S_20_03515","S_20_01019","S_20_00751","S_20_00948","S_20_02539")  # PAN

length(id.fst.5Q)
Q.res %>% subset(ID_GQ %in% id.fst.5Q) %>% nrow()
Q.res %>% subset(ID_GQ %in% id.fst.5Q) %>% group_by(Region1) %>% summarise(N = n())
Q.res %>% subset(ID_GQ %in% id.fst.5Q) %>% group_by(Membership.all) %>% summarise(N = n())  

gl.red.5Q <- gl.keep.ind(gl.final.5Q, id.fst.5Q, recalc = T, verbose = 5)  # keep only 35 individual of choice
table(pop(gl.red.5Q), useNA = "ifany")
FST.red.5Q <- gl.fst.pop(gl.red.5Q, nboots = 999, percent = 95, nclusters = 40)  # using 9999 results are pretty much the same


# save(list = c("FST.5Q","FST.red.5Q"),
#      file = file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_all_K5.Rdata"))
load(file.path("./02_Results/01_ddRAD_Bringloe/03_FST/Fst_all_K5.Rdata"))

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



## Figures: SLE vs CSB vs JAM-BEL vs LGR vs PAN ---------------------------

FST.table.5Q <- FST.5Q %>% table.fst()
FST.table.red.5Q <- FST.red.5Q %>% table.fst()

FST.5Q %>% table.fst() %>%
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst),
            Max.Pvalue = round(max(`p-value`),))
#       Mean        sd        Min       Max Max.Pvalue
# 0.05169948 0.0370715 0.01188907 0.1046736          0

FST.red.5Q %>% table.fst() %>%
  summarise(Mean = mean(Fst),
            sd = sd(Fst),
            Min = min(Fst),
            Max = max(Fst),
            Max.Pvalue = round(max(`p-value`),))
#       Mean         sd        Min       Max Max.Pvalue
# 0.06235071 0.04198957 0.01062719 0.1332582          0


FST.ddRAD.5Q <- FST.5Q %>% heat.fst() %>% 
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)
write.csv(FST.ddRAD.5Q, file = "./02_Results/01_ddRAD_Bringloe/03_FST/03a_all_SNPs/FST_ddRAD_all.csv", row.names = F)

FST.ddRAD.red.5Q <- FST.red.5Q %>% heat.fst() %>%
  select(Population1, Population2, Fst) %>% arrange(Population1, Population2)

gFST.ddRAD.5Q <-  FST.ddRAD.5Q %>%
  mutate(Population1 = factor(Population1, levels = c("PAN","JAM","LGR","CSB","SLE")),
         Population2 = factor(Population2, levels = c("PAN","JAM","LGR","CSB","SLE"))) %>% 
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
gFST.ddRAD.5Q




# Heterozygosity and Fis --------------------------------------------------

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity"))
}

if(!file.exists(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04a_all_SNPs"))){
  dir.create(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04a_all_SNPs"))
  print(file.path("./02_Results/01_ddRAD_Bringloe", "04_Heterozygosity", "04a_all_SNPs"))
}

vcf.path <- file.path("00_Data/03_ddRAD_Bringloe/populations.26019snps.638ind.final.recode.vcf")

cmd4 <- paste("--vcf", file.path(current.wd, vcf.path),
              "--het")
cmd4

setwd(file.path(current.wd, "02_Results/01_ddRAD_Bringloe/04_Heterozygosity/04a_all_SNPs/"))  # specify wd so that .ldepth.mean file is written in the right place

A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
tail(A4)

cat(file = "Heterozygosity.ddRAD.log",
    "\n", cmd4, "\n",
    A4, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

het.ddRAD <- read.delim("./02_Results/01_ddRAD_Bringloe/04_Heterozygosity/04a_all_SNPs/out.het", skip=0, sep = "\t", header = T )
het.ddRAD %>% head()

pop.het.ddRAD <- Q.res %>% left_join(het.ddRAD, by = c("ID_GQ"="INDV")) %>%
  select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Community, Membership.all, O.HOM., E.HOM., N_SITES, 'F') %>% 
  rename(Fis = 'F') %>% 
  mutate(Marker = "ddRAD")

pop.het.lcWGS <- read.csv("./02_Results/01_ddRAD_Bringloe/04_Heterozygosity/beluga_lcWGS_9v23.csv", stringsAsFactors = F)
pop.het.lcWGS %>% head()
pop.het.lcWGS <- pop.het.lcWGS %>% rename(Fis = 'F') %>% 
  mutate(Membership = str_replace(Membership, "CBS", "CSB"),
         Membership = str_replace(Membership, "JAM-BEL", "JAM")) %>%
  mutate(Marker = "lcWGS")

pop.het.ddRAD <- pop.het.ddRAD %>% mutate(ID = gsub("_rep", "", ID_GQ))
length(which(pop.het.ddRAD$ID %in% pop.het.lcWGS$ID))  # 73 samples in both dataset
nrow(pop.het.ddRAD) + nrow(pop.het.lcWGS) - length(which(pop.het.ddRAD$ID %in% pop.het.lcWGS$ID))  # 905

pop.het.all <- pop.het.ddRAD %>% select(ID_GQ, Membership.all, O.HOM., E.HOM., N_SITES, Fis, Marker) %>% rename("Membership" = "Membership.all") %>% 
  bind_rows(pop.het.lcWGS %>% select(ID_GQ=ID, Membership, O.HOM., E.HOM., N_SITES, Fis, Marker))  # ID_GQ=ID renames ID column in pop.het.lcWGS to ID_GQ, so the merging ends well

gFIS.all <- pop.het.all %>% subset(Membership %nin% "WHB") %>% 
  ggplot(aes(y = as.numeric(as.character(Fis)),
             x = factor(Membership, levels = c("RES","PAN","JAM","LGR","CSB","SLE")),
             col = factor(Membership, levels = c("RES","PAN","JAM","LGR","CSB","SLE")),
             shape = Marker)) +
  ggforce::geom_sina(alpha = 0.4, size = 8) +
  # scale_color_manual(name = "Genetic cluster", values = c("#c8ebf6ff","royalblue1","#e5b65fff","#f6835f","#49c687ff","#79c9e4ff"),
  #                    labels = c("Resolute passage","Panmictic","James Bay","Little-Great Whale Rivers","Cumberland Sound","St. Lawrence estuary")) +
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
# pdf(file = "02_Results/01_ddRAD_Bringloe/04_Heterozygosity/FIS.lcWGS.ddRAD.230717.pdf", width = 20, height = 10)
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


wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
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
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 62, p-value = 1.665e-14
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 1580.5, p-value = 1.647e-08
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "ddRAD"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "ddRAD"],
            alternative = "two.sided")  # U = 12974, p-value = 0.2983
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "ddRAD"],
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
            pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 1123, p-value = 0.4315
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "SLE" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 233, p-value = 1.544e-10
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "CSB" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 3177, p-value = 0.9012
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
            pop.het.all$Fis[pop.het.all$Membership %in% "LGR" & pop.het.all$Marker %in% "lcWGS"],
            alternative = "two.sided")  # U = 2598, p-value = 0.272
wilcox.test(pop.het.all$Fis[pop.het.all$Membership %in% "PAN" & pop.het.all$Marker %in% "lcWGS"],
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




# Misc stats --------------------------------------------------------------

## No reads (after demultiplexing) ----------------------------------------

reads <- read.csv("../MOBELS_PopStructure_ddRAD_MS_SNPdiscovery/02_Results/00_Stacks/03_Demultiplex/AllIndividuals_Nreads.csv", stringsAsFactors = F)

lapply(unique(reads$Run), function(i){
  d <- reads[reads$Run %in% i,]
  s <- sum(d$Total)
  print(paste(i, s))
})

reads$Total.single <- reads$Total/2
reads.red <- reads %>% filter(Filename %in% indNames(gl.final)) %>% 
  group_by(Filename) %>% summarise(Total = sum(Total)/2)  # divided by two fr the pair-end
summary(reads.red$Total)


## Stats on ddRAD loci depth ----------------------------------------------

depth <- read.csv("../MOBELS_PopStructure_ddRAD_MS_SNPdiscovery/02_Results/00_Stacks/05b_Stacks.ref/01_Bringloe_22xii22/01_popgen/BelugaID_unique_noMarralik_RefGen_Bringloe_NreadsNloci.csv", stringsAsFactors = F)

depth.red <- depth %>% filter(sample %in%  indNames(gl.final))  # use first gl.final in filtering script
summary(depth$mean_cov_ns)


# SIDE PROJECT
# ## Table nDNA CSB beluga --------------------------------------------------
# 
# admix.62ind.MAF10NA05.neutral.k2 <- Q.res.62ind.MAF10NA05.neutral.K2 %>% pivot_longer(cols = c("Q1","Q2"), names_to = "Group", values_to = "Q") %>% 
#   mutate(Group = factor(Group, levels = c("Q1","Q2")))  # admixture membership probability
# 
# cs <- pop.data %>% filter(ID_GQ %in% indNames(gl.final)) %>% 
#   filter(Pop %in% "CSB") %>% 
#   select(ID_GQ, Annee_echantillonnage, Mois_echantillonnage, Jour_echantillonnage, Pop, Lieu_echantillonnage, Latitude_echantillonnage_DD,
#          Longitude_echantillonnage_DD, Sexe_laboratoire, Haplotype, Private) %>% 
#   left_join(subset(admix.62ind.MAF10NA05.neutral.k2, Group %in% "Q2", select = c(ID_GQ, Q))) %>% 
#   arrange(Annee_echantillonnage, Mois_echantillonnage, Jour_echantillonnage, ID_GQ)
# colnames(cs) <- c("ID","Year","Month","Day","Sampling.region","Sampling.site","Latitude","Longitude","qPCR.Sex","Dloop.Haplotype","Private.Dloop","Admixture.CSBgroup.membership")
# 
# write.csv(cs, file.path("./00_Data/00_Dataset/", "CSB.ddRAD.metadata.27ind.csv"), row.names = F, quote = F)
# 
# 
# 
# ## Check what did you do this for
# 
# csb <- Q.res.54ind.MAF10NA05.neutral %>%
#   left_join(pop.data[,c("ID_GQ","Pop","Lieu_echantillonnage","Latitude_echantillonnage_DD","Longitude_echantillonnage_DD")]) %>%
#   filter(Pop %in% "CSB")
# colnames(csb)[c(6,7)] <- c("lat","lon")
# 
# csb <- csb %>% 
#   mutate(area = ifelse(lon %in% -67.44740 & lat %in% 66.56740, "Clearwater Fiord",
#                        ifelse(lon < -67.24679 & lon > -67.77979 & lat < 66.37291 & lat > 66.06069, "Kingalo",
#                               ifelse(lon < -65.35248 & lon > -67.17770 & lat < 66.36532 & lat > 65.89349, "North Stratum",
#                                      ifelse(lon > -63.33732, NA,
#                                             "West Stratum"))))) %>% 
#   subset(!is.na(area))
# 
# csb %>% group_by(area) %>% summarise(Min = )
# 

