# Info --------------------------------------------------------------------
# 
# Authors: Audrey Bourret, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice-Lamontagne 
# Date: 2022-07-05
# 
# Overview: Filtering SNPs (obtained from alignment to Trevor's reference genome) pipeline with STACKS version >2 population
# Filteris steps:
# 1. populations pipeline: -r (0.75) and MAF (0.01)
# 2. vcftools: min read depth - tried 10, 12, and 15 (+ max DP 26.656)
# 3. vcftools: missing data
# 4. vcftools: minium ID coverage 5X
# 5. vcftools
# http://vcftools.sourceforge.net/man_latest.html
# http://www.ddocent.com/filtering/


# Libraries ---------------------------------------------------------------

library(parallel) # to detect the number of cores
library(tidyverse)
library(readxl)

library(vcfR)
library(adegenet)
library(hierfstat)

library(here)

# Internal functions
for(i in 1:length(list.files("./01_Codes/Functions"))){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}
`%nin%` <- Negate(`%in%`)

# Paths - stacks and plink
stack2.55_path <- "/home/genyoda/Documents/Programs/stacks-2.55/" 
plink_path <- "/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606/" 



# Data --------------------------------------------------------------------

stacks.ref.path.popgen <- get.value("stacks.ref.path.bringloe.popgen")
filter.ref.path <-  file.path(get.value("filter.ref.path.bringloe"))

pop.data <- read_csv(file.path(get.value("info.path"),"Project_Infos.csv"))
pop.data 

# What is the sample size
nrow(pop.data[pop.data$Espece %nin% "Monodon" & pop.data$Region_echantillonnage %in% c("CSB","NWH","REB","NHB","SWH"),])  # 187 Delphinapterus - 8 Monodon
pop.data %>% pull(Cat_sample) %>% table()
pop.data %>% pull(Numero_unique_specimen) %>% table()
pop.data %>% pull(Region_echantillonnage) %>% table()
pop.data %>% pull(Espece) %>% table()

pop.data %>% filter(Cat_sample == "Sample") %>% 
  group_by(Region_echantillonnage) %>% 
  summarise(N = n()) #%>% write_csv("clipboard")

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

numCores




# Coverage ----------------------------------------------------------------

list.files(get.value("stacks.ref.log.bringloe.popgen"))
cov.data <- read_csv(file.path(get.value("stacks.ref.log.bringloe.popgen"), "BelugaID_unique_noMarralik_RefGen_Bringloe_NreadsNloci.csv"))

cov.data %>%
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = Espece)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  facet_wrap(~Region_echantillonnage, nrow = 2) +
  labs(x = "Mean coverage") +
  theme_bw() +
  theme(legend.position = "bottom")

cov.data %>%
  # subset(Espece %in% "Delphinapterus") %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = Region_echantillonnage, col = as.factor(Region_echantillonnage))) +
  geom_violin(col = "black") +
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = c(5,10), lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") +
  theme_bw() +
  theme(legend.position = "none")

cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Espece))) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 60000, lty = "dashed", col = "darkgray") +
  #  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")

cov.data %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.factor(Annee_echantillonnage), col = as.factor(Espece))) +
  geom_boxplot(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")


# List of individuals to use for next steps 
# Duplicates have already been removed
# For now, I'm keeping IDs with mean_cov_ns < 5
# ID.coverage <- cov.data %>%
#   filter(as.numeric(as.character(mean_cov_ns)) >= 5) %>%  # use mean_cov_ns (weighted mean coverage) to put more weight on loci shared among multiple samples
#   mutate(Pop = "NoPop") %>% 
#   select(sample, Pop)

cov.data.uni <- cov.data %>% 
  group_by(Numero_unique_specimen) %>%  # group by id to find duplicates
  slice_max(as.numeric(as.character(mean_cov_ns)))  %>%  # among duplicated specimens, remove those with lower reads coverage
  ungroup()  # to remove grouping variable Numero_unique_specimen: otherwise ID.coverage will contain it as column even when selecting for sample and Pop only

ID.coverage <- cov.data.uni %>%
  filter(as.numeric(as.character(mean_cov_ns)) >= 5, Espece %in% "Delphinapterus") %>% 
  mutate(Pop = "NoPop") %>% 
  select(sample, Pop)


cov.data %>% group_by(as.numeric(as.character(mean_cov_ns)) >= 5) %>% summarise(N = n())
cov.data %>% group_by(Region_echantillonnage) %>% summarise(N = n())
cov.data %>% filter(as.numeric(as.character(mean_cov_ns)) >= 5) %>% group_by(Region_echantillonnage) %>% summarise(N = n())

if(!file.exists(file.path(filter.ref.path))){
  dir.create(file.path(filter.ref.path))
  print(file.path(filter.ref.path))
}

write.table(ID.coverage, 
            file = file.path(filter.ref.path, "popmap.coverage.popgen5X.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)




# Filtering ---------------------------------------------------------------

# Parameters
r.value       <- 0.75  # Minimum within pop
R.value       <- 0.75  # Minimum overall
maf.value     <- 0.01  # Overall MAF
n.pop         <- ID.coverage %>% pull(Pop) %>% unique() %>% length()  # treated 
#maf.pop.value <- 0.05




## Filtering #1: -r and -MAF ----------------------------------------------

if(!file.exists(file.path(filter.ref.path, "01a_r75_MAF01"))){
  dir.create(file.path(filter.ref.path, "01a_r75_MAF01"))
  print(file.path(filter.ref.path, "01a_r75_MAF01"))
}

cmd <- paste("-P", stacks.ref.path.popgen, 
             "-M", file.path(filter.ref.path, "popmap.coverage.popgen5X.txt"),
             "--out-path", file.path(filter.ref.path, "01a_r75_MAF01"),
             "-t", 20,  # done with 20 on 230130
             "-r", r.value,           
             "-R", R.value, 
             "--min-maf", maf.value,
             "--min-populations", n.pop,  # treated like if all samples came from one population for filtering purposes (but see O'Leary et al. 2018, p. 3197)
             #"--smooth",
             #"--write-single-snp",
             "--vcf",
             "--plink"
)
cmd

A <- system2(paste0(stack2.55_path, "populations"), cmd, stdout=T, stderr=T)
A

# If you want the big file to be ignored, run the following :

cat("*.tsv", "*.vcf", "*.distribs", "*.map", "*.ped", "*.RData", ".gitignore", "*.log", sep = "\n",
    file = file.path(filter.ref.path, "01a_r75_MAF01", ".gitignore"))


## General check : Ho and Fis, He vs Ho outliers --------------------------
# Done using data post population command filtration (MAF01 + NA05 + without IDs < 5X)

# Load the post-filtration statistic table

filter.stat <- read.delim(file.path(filter.ref.path, "01a_r75_MAF01", "populations.sumstats.tsv"), 
                          skip=n.pop, sep = "\t", header = T )
names(filter.stat)[1] <- "Locus.ID"
nrow(filter.stat[unique(filter.stat$Locus.ID),])  # 88,544 unique loci after populations filtration

# nrow(filter.stat) / 15  # when specifying population - add 1 row per pop (15 was no of pop used by AB)

summary(filter.stat)
head(filter.stat)

# For fun, the distribution of the number of SNPs, by scaffold and locus

filter.stat %>% group_by(Locus.ID) %>% summarise(Nsnps = n()) %>% 
  group_by(Nsnps) %>% summarise(Nloc = n()) %>% View()

filter.stat %>% group_by(Chr) %>% summarise(Nloc = length(Locus.ID %>% unique())  ) %>% 
  group_by(Nloc) %>% summarise(Nscaffold = n()) %>% View()

# Check the distribution of Fis and Ho -> possibly of those with very high number of alleles detected

filter.stat %>% ggplot(aes(x = Obs.Het, fill=Pop.ID)) +
  geom_histogram() +
  theme_bw()

filter.stat %>% ggplot(aes(x = Fis)) +
  geom_histogram() +
  theme_bw()


filter.stat %>% ggplot(aes(x = Obs.Het, y = Fis, col = Pop.ID)) +
  #geom_point(alpha = 1/100) +
  geom_point()+
  geom_vline(xintercept = 0.6, col = "red", lty = "dashed")+
  geom_hline(yintercept = c(-0.3), col = "red", lty = "dashed") +
  theme_bw()

# filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het, col = Fis)) +
#   geom_point() +
#   geom_vline(xintercept = 0.6, col = "red", lty = "dashed")+
#   scale_colour_gradientn(colours=rainbow(4)) +
#   geom_abline(slope = 1 ) +
#   #geom_hline(yintercept = c(-0.3,0.3), col = "red", lty = "dashed") +
#   theme_bw()

graph0 <- filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het)) +
  geom_point(alpha = 0.075) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  geom_abline(slope = 2, col = "blue", lty = "dashed") +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  #labs(title = "Ho vs He overall")+ 
  labs(x = "Observed heterozygosity", y = "Expected heterozigosity") +
  theme_bw() +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
graph0

ggsave(filename = file.path(filter.ref.path,"01a_r75_MAF01", "He_Ho.136771snps.png"),
       plot = graph0,
       width = 10, height = 10, units = "in")



## Filtering #2: Loci DP --------------------------------------------------

# It can only be done on a vcf file, so you need to run stacks::populations first
# To consider: you might want to leave all IDs with low (<5X) coverage in previous step, as they will be removed anyway after
# Or maybe re-filter for IDs with low coverage after this step

# https://speciationgenomics.github.io/filtering_vcfs/
# A possibility is to use the --min-meanDP; now creating a list of SNPs to removeb based on min and max per locus meanDP
# Previously: since the distribution of mean (or even median) DP of loci are skewed, it is a better idea to use the
# medianDP which we can estimate manually in R (old codes in Filtering #5: Too much depth), since the median is more robust to represent the central tendency of 
# a skewed distribution. It could be interesting to compare the results though

if(!file.exists(file.path(filter.ref.path, "01b_DP"))){
  dir.create(file.path(filter.ref.path, "01b_DP"))
  print(file.path(filter.ref.path, "01b_DP"))
}

vcf.path <- file.path(filter.ref.path, "01a_r75_MAF01", "populations.snps.vcf")


# Mean coverage depth per sample

cmd1a <- paste("--vcf", file.path(current.wd, vcf.path), 
               "--depth"  # Generates a file containing the mean depth per individual. This file has the suffix ".idepth"
)
cmd1a

setwd(file.path(filter.ref.path, "01b_DP"))  # specify wd so that .ldepth.mean file is written in the right place

A1a <- system2("vcftools", cmd1a, stdout=T, stderr=T)
A1a

cat(file = "populations.filt_GLOBAL_ind-mean-depth.log",
    "\n", cmd1a, "\n",
    A1a, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

idepth <- read.delim(file.path(filter.ref.path, "01b_DP", "out.idepth"), skip=0, sep = "\t", header = T )
idepth %>% head()

idepth %>% filter(MEAN_DEPTH >= 5) %>% nrow()  # 680 samples with mean depth >= 5 (AFTER R and MAF filtration)
cov.data %>% filter(mean_cov_ns >= 5) %>% nrow()  # 680 samples with mean depth >= 5 (BEFORE R and MAF filtration)
median(idepth$MEAN_DEPTH)  # 16.1868
median(cov.data$mean_cov_ns)  # 12.942


# Mean coverage depth per site across all samples

cmd1b <- paste("--vcf", file.path(current.wd, vcf.path), 
               "--site-mean-depth"  # Generates a file containing the mean depth per site averaged across all individuals
)
cmd1b

setwd(file.path(filter.ref.path, "01b_DP"))  # specify wd so that .ldepth.mean file is written in the right place

A1b <- system2("vcftools", cmd1b, stdout=T, stderr=T)
A1b

cat(file = "populations.filt_GLOBAL_site-mean-depth.log",
    "\n", cmd1b, "\n",
    A1b, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

ldepth <- read.delim(file.path(filter.ref.path, "01b_DP", "out.ldepth.mean"), skip=0, sep = "\t", header = T )
ldepth %>% head()
# ldepth$MEAN_DEPTH %>% hist()

median(ldepth$MEAN_DEPTH)  # 17.5459
ldepth %>% ggplot(aes(MEAN_DEPTH)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  xlim(0,50) +
  geom_vline(xintercept = median(ldepth$MEAN_DEPTH), color = "red", size = 1.25, lty = "dashed") +  # using median depth (mean = 20.something)
  theme_light()

filter.stat %>% mutate(Rownumber = row_number()) %>% 
  left_join(ldepth %>% mutate(Rownumber = row_number())) %>%
  # mutate(ChrT = Chr == CHROM,  # manual check to see if Chromosome name and position are the same since there is no LocusID in ldepth
  #        PosT = BP == POS)
  # table(graph0.2$ChrT)
  # table(graph0.2$PosT)
  summarise(N = table(MEAN_DEPTH > 100))  # 907 SNPs with depth > 50 (407 > 100)

graph1.0 <- filter.stat %>% mutate(Rownumber = row_number()) %>% 
  left_join(ldepth %>% mutate(Rownumber = row_number())) %>%
  filter(MEAN_DEPTH <= 50) %>% 
  ggplot(aes(x = Locus.ID, y = Fis, color = MEAN_DEPTH)) +
  geom_point(alpha = 0.25) +
  viridis::scale_color_viridis(discrete = F) +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = mean(filter.stat$Fis), col = "black", lty = "dashed") +
  #labs(title = "Ho vs He overall")+
  #labs(x = "Observed heterozygosity", y = "Expected heterozigosity") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
graph1.0

ggsave(filename = file.path(filter.ref.path,"01b_DP", "Fis.MeanDepthPerSNP.png"),
       plot = graph1.0,
       width = 15.25, height = 10, units = "in")


# per-SNP quality ("--site-quality") doesn't work. I think it's because NAs are present in multiple 'gt_GQ' rows
# gstacks already had a default filter for reads quality (Q = 10)


vcf.data <- vcfR::read.vcfR(vcf.path)
gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # nrow = 99557430: multiple alleles for the same locus
# length(unique(gt.tidy$Key))  # 133455

vcf.fix <- as.data.frame(vcf.data@fix) %>% mutate(Key = 1:nrow(.))


# save(list = c("gt.tidy", "vcf.fix"),
#      file = file.path(filter.ref.path,"01b_DP", "SNPs_DP.RData"))
load(file.path(filter.ref.path, "01b_DP","SNPs_DP.RData"))


# group by snps
gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>%
  left_join(vcf.fix %>% select(Key, ID))

gt.key %>% ggplot(aes(meanDP)) +  # exactly the same distribution (phew) computed in LL 293-298
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  xlim(0,50) +
  theme_light()

graph1.0 <- gt.key %>% 
  ggplot(aes(x = medianDP, y = meanDP, col = maxDP)) +
  geom_jitter(alpha = 0.5)+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_color_viridis_c() +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graph1.0

ggsave(filename = file.path(filter.ref.path,"01b_DP", "DepthCoverage.png"),
       plot = graph1.0,
       width = 10, height = 10, units = "in")

gt.key %>% ggplot(aes(x = medianDP)) +
  geom_histogram()+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()

quantile(gt.key$meanDP, .99)  # 36.83506 (median = 20)
quantile(gt.key$meanDP, .50)  # 17.54586 (median = 13)
LOC.DP.MIN.10 <- gt.key %>% filter(meanDP <= 10) %>% pull(ID)  # 5671 loci
LOC.DP.MIN.12 <- gt.key %>% filter(meanDP <= 12) %>% pull(ID)  # 10192 loci - didn't seem to make much of a difference with LOC.DP.MIN10 when looking to He vs Ho plot
LOC.DP.MIN.15 <- gt.key %>% filter(meanDP <= 15) %>% pull(ID)  # 28646 loci
# LOC.DP.MIN <- gt.key %>% filter(medianDP <= 12) %>% pull(ID)  # 39357 loci (was 39372 loci using TB's ref geno; was 12533 with medianDP < 10 and 27925 with < 12)
# LOC.DP.MAX <- gt.key %>% filter(medianDP > 28) %>% pull(ID)  # 99th percentile - 1339 loci (1280 loci using 29 for Trevor's september ref genome)
# Increased to 12 to exclude loci with higher homozygosity than expected
# For ms/resdoc: we filterd loci by removing those with a median depth below 12 and higher than 29. The maximum value was chosen to exclude those loci that
# have too great a coverage meaning blah blah. For minim value, we tested a median depth estimated for each locus across all individuals of 10 and 15, but 
# proceeded with 12 because XX% of loci showed greater observed homozygosity than expected (after several filtering spets), which was highly likely caused by the low 
# coverage of certain alleles
LOC.DP.MAX <- gt.key %>% filter(meanDP > 36.83506) %>% pull(ID)  # 99th percentile - 1368 loci

LOC.ALL <- gt.key %>% pull(ID) %>% data.frame() %>% rename(ID='.')
LOC.DP.10 <- LOC.ALL %>% filter(ID %nin% c(LOC.DP.MAX,LOC.DP.MIN.10))  # 129732 loci
LOC.DP.12 <- LOC.ALL %>% filter(ID %nin% c(LOC.DP.MAX,LOC.DP.MIN.12))  # 125211 loci
LOC.DP.15 <- LOC.ALL %>% filter(ID %nin% c(LOC.DP.MAX,LOC.DP.MIN.15))  # 106757 loci

nrow(LOC.DP.10) + length(LOC.DP.MAX) + length(LOC.DP.MIN.10) == nrow(LOC.ALL)
nrow(LOC.DP.12) + length(LOC.DP.MAX) + length(LOC.DP.MIN.12) == nrow(LOC.ALL)
nrow(LOC.DP.15) + length(LOC.DP.MAX) + length(LOC.DP.MIN.15) == nrow(LOC.ALL)

write.csv(LOC.DP.MIN.10,file.path(filter.ref.path,"01b_DP", "Loc.MIN10.all.DP.csv"), row.names = F, quote = F)
write.csv(LOC.DP.MIN.12,file.path(filter.ref.path,"01b_DP", "Loc.MIN12.all.DP.csv"), row.names = F, quote = F)
write.csv(LOC.DP.MIN.15,file.path(filter.ref.path,"01b_DP", "Loc.MIN15.all.DP.csv"), row.names = F, quote = F)
write.csv(LOC.DP.MAX,file.path(filter.ref.path,"01b_DP", "Loc.MAX.all.DP.csv"), row.names = F, quote = F)

write.csv(LOC.DP.10,file.path(filter.ref.path,"01b_DP", "Loc.MIN10MAX.all.DP.csv"), row.names = F, quote = F)
write.csv(LOC.DP.12,file.path(filter.ref.path,"01b_DP", "Loc.MIN12MAX.all.DP.csv"), row.names = F, quote = F)
write.csv(LOC.DP.15,file.path(filter.ref.path,"01b_DP", "Loc.MIN15MAX.all.DP.csv"), row.names = F, quote = F)
LOC.DP.10 <- read.csv(file.path(filter.ref.path,"01b_DP", "Loc.MIN10MAX.all.DP.csv"), stringsAsFactors = F)
LOC.DP.12 <- read.csv(file.path(filter.ref.path,"01b_DP", "Loc.MIN12MAX.all.DP.csv"), stringsAsFactors = F)
LOC.DP.15 <- read.csv(file.path(filter.ref.path,"01b_DP", "Loc.MIN15MAX.all.DP.csv"), stringsAsFactors = F)


# Do the filtrations

cmd2a <- paste("--vcf", file.path(vcf.path),
               "--recode",
               #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
               #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
               "--snps", file.path(current.wd,filter.ref.path,"01b_DP", "Loc.MIN10MAX.all.DP.csv"),
               "--out", file.path(current.wd,filter.ref.path,"01b_DP", paste0("populations.",LOC.DP.10 %>% nrow(),"snps.",length(ID.coverage$sample),"ind.all"))
)
cmd2a  # populations.127118snps.746ind.all.recode.vcf

A2a <- system2("vcftools", cmd2a, stdout=T, stderr=T)
A2a

cat(file = file.path(current.wd,filter.ref.path,"01b_DP","VCFtools_SNPsMin10Max36.DP.log"),
    "\n", cmd2a, "\n",
    A2a, # what to put in my file
    append= F, sep = "\n")


cmd2b <- paste("--vcf", file.path(vcf.path),
               "--recode",
               #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
               #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
               "--snps", file.path(current.wd,filter.ref.path,"01b_DP", "Loc.MIN12MAX.all.DP.csv"),
               "--out", file.path(current.wd,filter.ref.path,"01b_DP", paste0("populations.",LOC.DP.12 %>% nrow(),"snps.",length(ID.coverage$sample),"ind.all"))
)
cmd2b  # populations.121831snps.746ind.all.recode.vcf

A2b <- system2("vcftools", cmd2b, stdout=T, stderr=T)
A2b

cat(file = file.path(current.wd,filter.ref.path,"01b_DP","VCFtools_SNPsMin12Max36.DP.log"),
    "\n", cmd2b, "\n",
    A2b, # what to put in my file
    append= F, sep = "\n")


cmd2c <- paste("--vcf", file.path(vcf.path),
               "--recode",
               #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
               #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
               "--snps", file.path(current.wd,filter.ref.path,"01b_DP", "Loc.MIN15MAX.all.DP.csv"),
               "--out", file.path(current.wd,filter.ref.path,"01b_DP", paste0("populations.",LOC.DP.15 %>% nrow(),"snps.",length(ID.coverage$sample),"ind.all"))
)
cmd2c  # populations.98059snps.746ind.all.recode.vcf

A2c <- system2("vcftools", cmd2c, stdout=T, stderr=T)
A2c

cat(file = file.path(current.wd,filter.ref.path,"01b_DP","VCFtools_SNPsMin15Max36.DP.log"),
    "\n", cmd2c, "\n",
    A2c, # what to put in my file
    append= F, sep = "\n")


# Verify differences between alternative minDP filters
list.files(file.path(filter.ref.path, "01b_DP"))

vcf.path10 <- file.path(filter.ref.path, "01b_DP", "populations.129732snps.680ind.all.recode.vcf")
vcf.path12 <- file.path(filter.ref.path, "01b_DP", "populations.125211snps.680ind.all.recode.vcf")
vcf.path15 <- file.path(filter.ref.path, "01b_DP", "populations.106757snps.680ind.all.recode.vcf")

l.path <- list(vcf.path10,vcf.path12,vcf.path15)
l.vcf.data <- lapply(l.path, function(i) vcfR::read.vcfR(i))  # import vcf files
l.gi.data <- lapply(l.vcf.data, function(i) vcfR::vcfR2genind(i))  # transfrom vcfs to genind objects
l.div <- lapply(l.gi.data, function(i) summary(i))  # summary of genind object

# save(list = c("l.div","l.vcf.data"),
#      file = file.path(filter.ref.path,"01b_DP", "minDP.RData"))
load(file.path(filter.ref.path, "01b_DP","minDP.RData"))

SCAFFOLD.info <- l.vcf.data[[3]]@fix %>% as.data.frame() %>%  # to create column RADloc and check how many RAD loci are kept after the filtering
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 69936

l.div.graph <- lapply(l.div, function(i) data.frame(ID = names(i$Hobs),
                                                    Hobs = i$Hobs,
                                                    Hexp = i$Hexp))

div10.graph <- l.div.graph[[1]]
div10.graph %>% ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.075) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  geom_abline(slope = 2, col = "blue", lty = "dashed") +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  #labs(title = "Ho vs He overall")+ 
  labs(x = "Observed heterozygosity", y = "Expected heterozigosity") +
  theme_bw() +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
div12.graph <- l.div.graph[[2]]
div15.graph <- l.div.graph[[3]]

l.graph.div <- lapply(l.div.graph, function(g){
  gg <- g %>% ggplot(aes(x = Hobs, y = Hexp)) +
    geom_point(alpha = 0.075) +
    scale_colour_distiller(palette = "Spectral") +
    geom_vline(xintercept = 0.6) +
    geom_abline(slope = 2, col = "blue", lty = "dashed") +
    geom_abline(slope = 1, col = "red", lty = "dashed") +
    labs(x = "", y = "") +
    theme_bw() +
    theme(axis.text = element_text(size = 15, colour = "black"),
          axis.title = element_text(size = 18, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  gg
})

graph.minDP10 <- l.graph.div[[1]]
graph.minDP10
graph.minDP12 <- l.graph.div[[2]]
graph.minDP12
graph.minDP15 <- l.graph.div[[3]]
graph.minDP15

graph.minDP <- filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het)) +
  geom_point(alpha = 0.075) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  geom_abline(slope = 2, col = "blue", lty = "dashed") +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


graph.minDP.all <- cowplot::plot_grid(graph.minDP, graph.minDP10, graph.minDP12, graph.minDP15,
                                      nrow = 2, ncol = 2, labels = "AUTO", label_x = 0.075, label_y = 0.985) +
  cowplot::draw_label("Expected heterozygosity", x = 0, y = 0.5, vjust = 1.5, angle = 90, size = 18) +
  cowplot::draw_label("Observed heterozygosity", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = 18)
graph.minDP.all

# minDP = 15 seems to be the one working best (kept 78.05% of 'original' SNPs)

ggsave(filename = file.path(filter.ref.path, "01b_DP", "He_Ho_minDP_20230418.png"),
       plot = graph.minDP.all,
       width = 20, height = 20, units = "in")

cat("*.vcf", "*.RData", "!.gitignore", sep = "\n",
    file = file.path(filter.ref.path, "01b_DP", ".gitignore"))


# Verify SNP depth for populations.106757snps.680ind.all.recode.vcf

vcf.data <- l.vcf.data[[3]]

gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))  # nrow = 72594760: multiple alleles for the same locus
# length(unique(gt.tidy$Key))  # 106757

vcf.fix <- as.data.frame(vcf.data@fix) %>% mutate(Key = 1:nrow(.))

# group by snps
gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>%
  left_join(vcf.fix %>% select(Key, ID))

gt.key %>% ggplot(aes(meanDP)) +  # although I filtered for meadianDP < 29, some meanDP are still quite high
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

graph1.1 <- gt.key %>% 
  ggplot(aes(x = medianDP, y = meanDP, col = maxDP)) +
  geom_jitter(alpha = 0.5)+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_color_viridis_c() +
  labs(title = "SNPs coverage depth")+
  theme_bw()
graph1.1

ggsave(filename = file.path(filter.ref.path,"01b_DP", "DepthCoverage1.png"),
       plot = graph1.1 ,
       width = 10, height = 10, units = "in")

gt.key %>% ggplot(aes(x = medianDP)) +
  geom_histogram()+
  geom_vline(xintercept = quantile(gt.key$medianDP, .99), lty = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(title = "SNPs coverage depth")+
  theme_bw()

cat("*.vcf", "*.RData", "!.gitignore", sep = "\n",
    file = file.path(filter.ref.path, "01b_DP", ".gitignore"))



## Filtration #3: Missing data --------------------------------------------

# Verification: loci with high proportion of missing alleles

list.files(file.path(filter.ref.path, "01b_DP"))

vcf.path <- file.path(filter.ref.path, "01b_DP", "populations.106757snps.680ind.all.recode.vcf")

cmd3a <- paste("--vcf", file.path(current.wd, vcf.path), 
               #"--out",  "MAX2.NArm",
               #" --max-alleles", 2,
               #"--max-missing", 0.70,  # exclude sites on the base of the proportion of missing data (between 0 and 1 where 0 sites completely missing allowed; 
               # 1 no missing data allowed)
               "--missing-site"  # missingness on a per-site basis
               #"--missing-indv"
               #"--kept-sites"
               #"--recode"
)
cmd3a

if(!file.exists(file.path(filter.ref.path, "01c_MissingData"))){
  dir.create(file.path(filter.ref.path, "01c_MissingData"))
  print(file.path(filter.ref.path, "01c_MissingData"))
}

setwd(file.path(filter.ref.path, "01c_MissingData")) 

A3a <- system2("vcftools", cmd3a, stdout=T, stderr=T)
A3a

cat(file = "populations.filt_GLOBAL_5x_Loci_wMissing.log",
    "\n", cmd3a, "\n",
    A3a, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

lmiss <- read.delim(file.path(filter.ref.path, "01c_MissingData", "out.lmiss"), skip=0, sep = "\t", header = T )
lmiss %>% head()

graph2.0 <- lmiss %>% ggplot(aes(x=F_MISS)) + 
  geom_histogram() + 
  geom_vline(xintercept = 0.10, lty = "dashed", col = "red") +
  theme_bw()
graph2.0  

ggsave(filename = file.path(filter.ref.path, "01c_MissingData", "NAbyLocus_20230418.png"),
       plot = graph2.0,
       width = 4, height = 4, units = "in")

good.LOC <- lmiss %>% filter(F_MISS < .10) %>% select(CHR,POS)

nrow(good.LOC)  # 92,115 SNPs

# Do the filtration: loci with missing proportion > 10%

cmd3b <- paste("--vcf", file.path(current.wd, vcf.path), 
               "--recode",
               "--max-missing", "0.9",  # Exclude sites on the basis of the proportion of missing data (0 allows completely missing sites, 1, no missing data allowed)
               #"--snps", file.path(current.wd,get.value("filter.ref.path"),"A_GLOBAL_PB_5x", "01c_Indiduals_wMissing","LOC.10.csv"),  # I could have used good.LOC to filter snps
               "--out", vcf.path %>% str_replace("01b_DP", "01c_MissingData") %>%
                 str_replace("populations.106757snps.680ind.all.recode.vcf", paste0("populations.",nrow(good.LOC),"snps.","680ind.all"))
)
cmd3b  # populations.74681snps.746ind.all (was populations.87600snps.691ind.all with previous analyses with Jones)

A3b <- system2("vcftools", cmd3b, stdout=T, stderr=T)
A3b

cat(file = file.path(filter.ref.path, "01c_MissingData","populations.filt_GLOBAL_5x_b_Remove_Loci.log"),
    "\n", cmd3b, "\n",
    A3b, # what to put in my file
    append= F, sep = "\n")

# If you want the big file to be ignored, run the following :
cat("*.vcf", "!.gitignore", sep = "\n",
    file = file.path(filter.ref.path, "01c_MissingData", ".gitignore"))


# Verification: individual with > 30% missing data

list.files(file.path(filter.ref.path, "01c_MissingData"))

cmd3c <- paste("--vcf", file.path(current.wd, filter.ref.path, "01c_MissingData", "populations.92115snps.680ind.all.recode.vcf"), 
               #"--out",  "MAX2.NArm",
               #" --max-alleles", 2,
               #"--max-missing", 0.80,
               #"--missing-site",
               "--missing-indv"
               #"--kept-sites"
               #"--recode"
)
cmd3c

setwd(file.path(filter.ref.path, "01c_MissingData")) 

A3c <- system2("vcftools", cmd3c, stdout=T, stderr=T)
A3c

cat(file = "populations.filt_GLOBAL_5x_Individuals_wMissing.log",
    "\n", cmd3c, "\n",
    A3c, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

# One ind with more 30% missing value

imiss <- read.delim(file.path(filter.ref.path, "01c_MissingData", "out.imiss"), skip=0, sep = "\t", header = T )
imiss %>% head()
imiss %>% filter(F_MISS >.30)

imiss %>% left_join(pop.data %>% select(INDV = ID_GQ, Espece)) %>%  ggplot(aes(x=F_MISS, fill = Espece)) + geom_histogram()

graph2.1 <- imiss %>% left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
  ggplot(aes(x = Region_echantillonnage, y = F_MISS)) + 
  geom_boxplot(col = "black", fill = NA) +
  geom_jitter(height = 0, alpha = 0.2, col = "blue") +
  geom_hline(yintercept = 0.30, lty = "dashed", col = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
graph2.1

ggsave(filename = file.path(filter.ref.path, "01c_MissingData", "NAbyRegion_20230419.png"),
       plot = graph2.1,
       width = 6, height = 4, units = "in")


good.ID <- imiss %>% filter(F_MISS <= .30) %>% pull(INDV) %>% as.character()

pop.data %>% filter(ID_GQ %in% good.ID) %>% #select(ID_GQ, Gen_ZONE) %>% 
  group_by(Region_echantillonnage) %>% summarise(N = n())

#write.table(pop.info %>% filter(Sample %in% good.ID) 
#             %>% select(Sample, POP), 
#             file = file.path(get.value("info.path"), "popmap.imiss.txt"),
#             quote = FALSE, sep = "\t",
#             row.names = F, col.names = F)

cmd3d <- paste("--vcf", file.path(current.wd,filter.ref.path,"01c_MissingData","populations.92115snps.680ind.all.recode.vcf"), 
               "--recode",
               paste("--indv", good.ID, collapse = " "),
               "--out", file.path(current.wd, filter.ref.path, "01c_MissingData", paste0("populations.92115snps.",length(good.ID),"ind.all")
               )
)
cmd3d

A3d <- system2("vcftools", cmd3d, stdout=T, stderr=T)
A3d %>% tail()

cat(file = file.path(filter.ref.path, "01c_MissingData","VCFtools_RemoveIndMissing.log"),
    "\n", cmd3d, "\n",
    A3d, # what to put in my file
    append= F, sep = "\n")

## Filtration #4: Sample DP -----------------------------------------------

# Mean coverage depth per sample

list.files(file.path(filter.ref.path, "01c_MissingData"))

vcf.path <- file.path(filter.ref.path, "01c_MissingData", "populations.92115snps.679ind.all.recode.vcf")

cmd4a <- paste("--vcf", file.path(current.wd, vcf.path), 
               "--depth"  # Generates a file containing the mean depth per individual. This file has the suffix ".idepth"
)
cmd4a

if(!file.exists(file.path(filter.ref.path, "01d_DPsamples"))){
  dir.create(file.path(filter.ref.path, "01d_DPsamples"))
  print(file.path(filter.ref.path, "01d_DPsamples"))
}

setwd(file.path(filter.ref.path, "01d_DPsamples"))  # specify wd so that .idepth file is written in the right place

A4a <- system2("vcftools", cmd4a, stdout=T, stderr=T)
A4a

cat(file = "populations.filt_GLOBAL_postmissing_ind-mean-depth.log",
    "\n", cmd4a, "\n",
    A4a, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

idepth <- read.delim(file.path(filter.ref.path, "01d_DPsamples", "out.idepth"), skip=0, sep = "\t", header = T )
idepth %>% head()

idepth %>% filter(MEAN_DEPTH >= 5) %>% nrow()  # 677 samples with mean depth >= 5
cov.data %>% filter(mean_cov_ns >= 5) %>% nrow()  # 680 samples with mean depth >= 5 (BEFORE R and MAF filtration)
median(idepth$MEAN_DEPTH)  # 15.532
median(cov.data$mean_cov_ns)  # 12.942

idepth %>% left_join(pop.data %>% select(INDV = ID_GQ, Espece)) %>%  ggplot(aes(x=MEAN_DEPTH, fill = Espece)) + geom_histogram()

graph4.1 <- idepth %>% left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
  ggplot(aes(x = Region_echantillonnage, y = MEAN_DEPTH)) + 
  geom_boxplot(col = "black", fill = NA) +
  geom_jitter(height = 0, alpha = 0.2, col = "blue") +
  geom_hline(yintercept = 0.30, lty = "dashed", col = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
graph4.1

ggsave(filename = file.path(filter.ref.path, "01d_DPsamples", "SampleDPbyRegion_20230419.png"),
       plot = graph4.1,
       width = 6, height = 4, units = "in")


# Filtration: remove samples with coverage < 5X

good.ID <- idepth %>% filter(MEAN_DEPTH > 5) %>% pull(INDV) %>% as.character()

pop.data %>% filter(ID_GQ %in% good.ID) %>% #select(ID_GQ, Gen_ZONE) %>% 
  group_by(Region_echantillonnage) %>% summarise(N = n())

cmd4b <- paste("--vcf", file.path(current.wd,filter.ref.path,"01c_MissingData","populations.92115snps.679ind.all.recode.vcf"), 
               "--recode",
               paste("--indv", good.ID, collapse = " "),
               "--out", file.path(current.wd, filter.ref.path, "01d_DPsamples", paste0("populations.92115snps.",length(good.ID),"ind.5X")
               )
)
cmd4b

A4b <- system2("vcftools", cmd4b, stdout=T, stderr=T)
A4b %>% tail()

cat(file = file.path(filter.ref.path, "01d_DPsamples","VCFtools_CoverageSampleDP5X.log"),
    "\n", cmd4b, "\n",
    A4b, # what to put in my file
    append= F, sep = "\n")

if(!file.exists(file.path(filter.ref.path, "01d_DPsamples", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01d_DPsamples", ".gitignore")) 
}



## Filtration #5: HW outliers ---------------------------------------------

if(!file.exists(file.path(filter.ref.path, "01e_HW"))){
  dir.create(file.path(filter.ref.path, "01e_HW"))
  print(file.path(filter.ref.path, "01e_HW"))
}

# This part can be long ... (data convertion to hierfstat)

list.files(file.path(filter.ref.path, "01d_DPsamples"))

# Load the VCF file

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path, "01d_DPsamples", "populations.92115snps.677ind.5X.recode.vcf"))

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after previous filtering step
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 62195

head(vcf.data)

vcf_field_names(vcf.data , tag = "FORMAT")

# Extract raw info from vcf file
gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))
head(gt.tidy)

gt.meta <- gt.tidy %>% group_by(Indiv) %>%  summarise(Nsnps = length(gt_GT[!is.na(gt_GT)]),
                                                      DP = mean(gt_DP, na.rm=T)) %>%
  left_join(pop.data, by = c("Indiv" = "ID_GQ"))

head(gt.meta)

gt.meta %>% ggplot(aes(x = Region_echantillonnage, y = Nsnps)) +
  #geom_violin(fill="#C0C0C0", adjust = 1, scale = "count", trim = T) +
  # geom_jitter(height = 0, alpha = 1/5) +
  geom_jitter(height = 0, alpha = 1/5) +
  geom_boxplot(alpha = 0) +
  #  stat_summary(fun.data=mean_sdl, geom = "pointrange", color = "black")+
  #scale_y_continuous(trans = scales::log2_trans(), breaks = c(1,10,100,1000, 10000))+
  labs(y = "N snps", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

gt.meta %>%   ggplot(aes(x = DP, y = Nsnps, col = No_plaque_envoi)) +
  geom_point()+
  #  stat_summary(fun.data=mean_sdl, geom = "pointrange", color = "black")+
  #scale_y_continuous(trans = scales::log2_trans(), breaks = c(1,10,100,1000, 10000))+
  labs(y = "N snps", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Conversion to genlight and genind R formats

gl.data  <- vcfR::vcfR2genlight(vcf.data)
gi.data  <- vcfR::vcfR2genind(vcf.data) 


# Computing diversity (Ho, He)

div <- summary(gi.data)
str(div)
div.graph <- data.frame(ID = names(div$Hobs),
                        Hobs = div$Hobs,
                        Hexp = div$Hexp)

div.graph %>% head()

# save(list = c("div", "div.graph"),
#     file = file.path(filter.ref.path,"01e_HW", "SNPs_5x_HW.RES.RData"))
load(file.path(filter.ref.path, "01e_HW","SNPs_5x_HW.RES.RData"))

graph5.0 <-  div.graph %>% ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.1) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  labs(title= "Ho vs He overall") + 
  theme_bw()
graph5.0

# ggsave(filename = file.path(filter.ref.path, "01e_HW", "He_Ho_20230419.png"),
#        plot = graph5.0,
#        width = 4, height = 4, units = "in")

div.graph %>% ggplot(aes(x = Hobs, y = Hexp)) +  # density plot for Ho vs He
  geom_density_2d_filled() +
  geom_abline(slope = 1, col = "black", lty = "dashed") 


graph5.1 <-  div.graph %>% ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.1) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  geom_abline(slope = 2, col = "blue", lty = "dashed") +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  labs(x = "", y = "")+ 
  theme_bw() +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
graph5.2 <- cowplot::plot_grid(graph.minDP, graph.minDP15, graph5.1,
                               nrow = 1, ncol = 3, labels = "AUTO", label_x = 0.075, label_y = 0.985) +
  cowplot::draw_label("Expected heterozygosity", x = 0, y = 0.5, vjust = 1.5, angle = 90, size = 18) +
  cowplot::draw_label("Observed heterozygosity", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = 18)
graph5.2

ggsave(filename = file.path(filter.ref.path,"01e_HW", "He_Ho_20230419.png"),
       plot = graph5.2,
       width = 25, height = 10, units = "in")


gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>%
  left_join(vcf.fix %>% select(Key, ID))


# Add info on loci

div.graph <- div.graph %>% mutate(Loc = sapply(str_split(ID, ":"), `[`, 1)) 

Nspns.tab <- div.graph %>% group_by(Loc) %>% 
  summarise(N = n())              

# Keep 1 snp / rad loci
list.h06.unique <- div.graph %>% filter(Hobs <=0.6)  %>% distinct(Loc, .keep_all = T) 

# Keep ALL
list.h06.all  <- div.graph %>% filter(Hobs <=0.6) 

nrow(list.h06.unique)  # 62069
nrow(list.h06.all)  # 91967

write.csv(list.h06.unique %>% select(ID), file.path(filter.ref.path,"01e_HW", "Loc.h06.unique.csv"), 
          row.names = F, quote = F)

write.csv(list.h06.all %>% select(ID), file.path(filter.ref.path,"01e_HW", "Loc.h06.all.csv"), 
          row.names = F, quote = F)


# CREATE VCF WITH UNIQUE
list.files(file.path(filter.ref.path, "01d_DPsamples"))

cmd5a <- paste("--vcf", file.path(filter.ref.path, "01d_DPsamples", "populations.92115snps.677ind.5X.recode.vcf"), 
               "--recode",
               "--snps", file.path(here::here(),filter.ref.path,"01e_HW", "Loc.h06.unique.csv"),
               "--out", file.path(here::here(),filter.ref.path,"01e_HW", paste0("populations.",list.h06.unique %>% nrow(),"snps.677ind.H06.single"))
)
cmd5a

A5a <- system2("vcftools", cmd5a, stdout=T, stderr=T)  # populations.58229snps.682ind.H06.single: only two alleles per locus
A5a

cat(file = file.path(here::here(),filter.ref.path,"01e_HW","VCFtools_SnpWhiteList.H06.log"),
    "\n", cmd5a, "\n",
    A5a, # what to put in my file
    append= F, sep = "\n")


cmd5b <- paste("--vcf", file.path(filter.ref.path, "01d_DPsamples", "populations.92115snps.677ind.5X.recode.vcf"), 
               "--recode",
               "--snps", file.path(here::here(),filter.ref.path,"01e_HW", "Loc.h06.all.csv"),
               "--out", file.path(here::here(),filter.ref.path,"01e_HW", paste0("populations.",list.h06.all %>% nrow(),"snps.677ind.H06.all"))
)
cmd5b

A5b <- system2("vcftools", cmd5b, stdout=T, stderr=T)    # populations.86343snps.682ind.H06.all: all alleles available per locus
A5b

cat(file = file.path(here::here(),filter.ref.path,"01e_HW","VCFtools_SnpWhiteList.H06.log"),
    "\n", cmd5b, "\n",
    A5b, # what to put in my file
    append= T, sep = "\n")


list.files(file.path(filter.ref.path, "01e_HW"))

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path, "01e_HW", "populations.91967snps.677ind.H06.all.recode.vcf"))
SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing SNPs with He > 0.6
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, "_"), `[`,2), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 62069


# Add a gitignore if necessary

if(!file.exists(file.path(filter.ref.path, "01e_HW", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01e_HW", ".gitignore")) 
}


# rapid PCA
file.path(current.wd,filter.ref.path,"01e_HW") %>% list.files()

vcf.path <- file.path(filter.ref.path,"01e_HW", "populations.62069snps.677ind.H06.single.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# subset

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}


# Function to create a list of loci, from a genind object

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
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

LOC.MAF10.NA01 <- filter.MAF.NA(gi.data, MAF.trs = 0.10, NA.trs = 0.01)
#LOC.MAF05.NA05 <- filter.MAF.NA(gi.data, MAF.trs = 0.05, NA.trs = 0.05)


pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Region_echantillonnage)

table(pop(gl.data))

pop(gi.data) <- data.frame(ID_GQ = indNames(gi.data)) %>% 
  left_join(pop.data) %>% pull(Region_echantillonnage)

table(pop(gi.data))


# Function to to count proportion of NA loci per individual from a genlight object

count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}


# save(list = c("gl.data", "gi.data","LOC.MAF05.NA05"),
#      file = here(filter.ref.path, "01e_HW", "gi.gl.populations.H06.single.Rdata"))
load(file.path(filter.ref.path, "01e_HW","gi.gl.populations.H06.single.Rdata"))


na.info <- data.frame(ID_GQ = indNames(gl.data),
                      NNA = count.ind.na.gl(gl.data[, locNames(gl.data) %in% LOC.MAF05.NA05]))

pca.test <- glPca(gl.data[, locNames(gl.data) %in% LOC.MAF05.NA05], center = TRUE, scale = FALSE,  
                  parallel = TRUE, n.core = 8, nf = 1000)

# save(list = c("pca.test", "na.info"),
#      file = here(filter.ref.path, "01e_HW", "PCA.Rdata"))
load(file.path(filter.ref.path, "01e_HW","PCA.Rdata"))


gPCA12 <- pca.test %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA12

ggsave(filename = file.path(filter.ref.path, "01e_HW", "PCA_PC1vs2_20230419.png"),
       plot = gPCA12,
       width = 10, height = 8, units = "in")

gPCA23 <- pca.test %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC2, y = score.PC3, col = Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)"),
    y = paste0("PC3 (", QuickPop::pca_var(pca.test)$p.eig[3] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA23

ggsave(filename = file.path(filter.ref.path, "01e_HW", "PCA_PC2vs3_20230419.png"),
       plot = gPCA23,
       width = 10, height = 8, units = "in")

gPCA34 <- pca.test %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4, col = Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC3 (", QuickPop::pca_var(pca.test)$p.eig[3] %>% round(3) *100, "%)"),
    y = paste0("PC4 (", QuickPop::pca_var(pca.test)$p.eig[4] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA34

ggsave(filename = file.path(filter.ref.path, "01e_HW", "PCA_PC3vs4_20230419.png"),
       plot = gPCA34,
       width = 10, height = 8, units = "in")


gPCA.NA <- pca.test %>% QuickPop::pca_scoretable(naxe = 4) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC2, y = score.PC3, col = NNA)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  scale_color_distiller(palette = "Spectral") +
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)"),
    y = paste0("PC3 (", QuickPop::pca_var(pca.test)$p.eig[3] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA.NA

ggsave(filename = file.path(filter.ref.path, "01e_HW", "PCA_PC2vs3_NA_20240419.png"),
       plot = gPCA.NA,
       width = 10, height = 8, units = "in")


# Ho vs He AFTER filtering for MAF and NA - suspicious loci

div.MAF05.NA05 <- summary(gi.data[loc = LOC.MAF05.NA05])

str(div.MAF05.NA05)

div.graph.MAF05.NA05 <- data.frame(ID = names(div.MAF05.NA05$Hobs),
                                   Hobs = div.MAF05.NA05$Hobs,
                                   Hexp = div.MAF05.NA05$Hexp)

div.graph.MAF05.NA05 %>% head()

graph4.3 <-  div.graph.MAF05.NA05 %>% ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.1) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  geom_vline(xintercept = c(0.08, 0.15 ,0.22), color = "blue") +
  geom_hline(yintercept = c(0.18, 0.29), color = "blue") +
  labs(title = "Ho vs He overall")+ 
  theme_bw()
graph4.3

ggsave(filename = file.path(filter.ref.path, "01e_HW", "He_Ho_MAF05_NA05_20230419.png"),
       plot = graph4.3,
       width = 6, height = 6, units = "in")

# Fis by SNP position
fis.MAF05.NA05.graph <- data.frame(ID = names(div.MAF05.NA05$Hobs),
                                   Hobs = div.MAF05.NA05$Hobs,
                                   Hexp = div.MAF05.NA05$Hexp)

gi.MAF05.NA05.single <- gi.data[loc = LOC.MAF05.NA05]
pop(gi.MAF05.NA05.single) <- rep("HB", nrow(gi.MAF05.NA05.single@tab))  # genind2hierfstat requires that gi.data.single has pop specified for all IDs (for now all IDs same pop)
hfstat.MAF05.NA05.single <- genind2hierfstat(gi.MAF05.NA05.single)
fis.hfstat.MAF05.NA05.single <- basic.stats(hfstat.MAF05.NA05.single, diploid = T, digits = 3)
# names(fis.hfstat.MAF05.NA05.single)
# fis.hfstat.MAF05.NA05.single$pop.freq[1:10]
fis.hfstat.MAF05.NA05.single$n.ind.samp[1:10]
fis.MAF05.NA05.graph <- fis.MAF05.NA05.graph %>% cbind(fis.hfstat.MAF05.NA05.single$perloc$Fis)
colnames(fis.MAF05.NA05.graph)[4] <- "Fis"

graph4.4 <- fis.MAF05.NA05.graph %>% ggplot(aes(x = ID, y = Fis)) +
  geom_point(alpha = 0.075, colour = "cornflowerblue") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = mean(filter.stat$Fis), col = "black", lty = "dashed") +
  #labs(title = "Ho vs He overall")+ 
  #labs(x = "Observed heterozygosity", y = "Expected heterozigosity") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
graph4.4

ggsave(filename = file.path(filter.ref.path,"01e_HW", "Fis.MAF05NA05.single.png"),  # radloci different order than Fis.MeanDepthPerSNP.139317snps.png (graph0.1)
       plot = graph4.4,
       width = 15.25, height = 10, units = "in")




## Filtration #6: Plate effect --------------------------------------------

file.path(current.wd,filter.ref.path,"01e_HW") %>% list.files()

vcf.path <- file.path(filter.ref.path,"01e_HW", "populations.91967snps.677ind.H06.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

gl.data  <- vcfR::vcfR2genlight(vcf.data) 

# RDA on plate

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% left_join(pop.data) %>% pull(No_plaque_envoi)

pop(gl.data) %>% table()

library(vegan)

## Full model
RDA.plate <- vegan::rda(tab(gl.data) ~ pop(gl.data))
# RDA.plate.corrected <- vegan::rda(tab(gl.data[, locNames(gl.data) %nin% outliers_plate$Loci]) ~ pop(gl.data))  # Test by AB to check if it worked

# tab(gl.data.sex)[1:5,1:5]

anova(RDA.plate)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = tab(gl.data.sex) ~ pop(gl.data.sex))
#           Df Variance      F Pr(>F)    
# Model      1     16.1 1.3092  0.001 ***
# Residual 647   7937.6                  
RsquareAdj(RDA.plate)$adj.r.squared  # 0.006315175

res.plate <- as.data.frame(scores(RDA.plate, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.plate, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

res.plate.corrected <- as.data.frame(scores(RDA.plate.corrected, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.plate.corrected, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

if(!file.exists(file.path(filter.ref.path, "01f_Plate"))){
  dir.create(file.path(filter.ref.path, "01f_Plate"))
  print(file.path(filter.ref.path, "01f_Plate"))
}

# save(list = c("RDA.plate", "res.plate"),
#      file = file.path(filter.ref.path,"01f_Plate", "SNPs_5x_Plate.RData"))
load(file.path(filter.ref.path, "01f_Plate","SNPs_5x_Plate.RData"))

screeplot(RDA.plate)

gg.plate.rda <- res.plate %>% mutate(Region_test = ifelse(Region_echantillonnage %in% c("JAM", "SLE", "SAN"), Region_echantillonnage, "Other")) %>% 
  ggplot(aes(x = RDA1, y = RDA2, col = No_plaque_envoi)) +
  facet_wrap(~Region_test) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point() +
  ggtitle("Before correction") +
    #geom_point(data = res2, colour = "grey70", cex = 1) +
  #geom_histogram() +
  #geom_density()+
  labs(x=c("RDA 1")) +  theme_bw()  


gg.plate.corrected <- res.plate.corrected %>% mutate(Region_test = ifelse(Region_echantillonnage %in% c("JAM", "SLE", "SAN"), Region_echantillonnage, "Other")) %>% 
  ggplot(aes(x = RDA1, y = RDA2, col = No_plaque_envoi)) +
  facet_wrap(~Region_test) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point() +
  #geom_point(data = res2, colour = "grey70", cex = 1) +
  #geom_histogram() +
  #geom_density()+
  ggtitle("After correction") +
  labs(x=c("RDA 1")) +  theme_bw()  

gg.plate.combi <- ggpubr::ggarrange(gg.plate.rda, gg.plate.corrected, common.legend = T)

ggsave(filename = file.path(filter.ref.path,"01f_Plate", "Plate_RDA_20230419.png"),  # radloci different order than Fis.MeanDepthPerSNP.139317snps.png (graph0.1)
       plot = gg.plate.combi,
       width = 9, height = 7, units = "in", bg = "white")

# 
# hist(RDA.plate$CCA$v[,1])
# 
# sd(abs(RDA.plate$CCA$v[,1]))*1.96
# 
# quantile(abs(RDA.plate$CCA$v[,1]), .99)
# 
# # https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
#   x[x < lims[1] | x > lims[2]]               # locus names in these tails
# }
# 
# load.rda <- scores(RDA.plate, choices=c(1:3), display="species")
# 
# hist(load.rda[,1], main="Loadings on RDA1")
# hist(load.rda[,2], main="Loadings on RDA2")
# hist(load.rda[,3], main="Loadings on RDA3") 
# 
# cand1 <- outliers(load.rda[,1],3) # 38
# cand2 <- outliers(load.rda[,2],3) # 69
# cand3 <- outliers(load.rda[,3],3) # 34
# 
# ncand <- length(cand1) + length(cand2) + length(cand3)
# ncand
# 
# plot(RDA.plate, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
# points(RDA.plate, display="species", pch=21, cex=1, col="gray32", scaling=3)
# #points(RDA.plate, display="species", pch=21, cex=1, scaling=3)
# text(RDA.plate, scaling=3, display="bp", col="#0868ac", cex=1)
# 
# 3221 / 91000


# Source Capblancq

rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}


library(robust)
library(qvalue)

rdadapt_plate <- rdadapt(RDA.plate, 3)

## P-values threshold after Bonferroni correction
thres_plate <- 0.05/length(rdadapt_plate$p.values)

## Identifying the loci that are below the p-value threshold
outliers_plate <- data.frame(Loci = colnames(tab(gl.data))[which(rdadapt_plate$p.values<thres_plate)], 
                                   p.value = rdadapt_plate$p.values[which(rdadapt_plate$p.values<thres_plate)])
outliers_plate    
outliers_plate  %>% nrow()

## List of outlier names

## Formatting table for ggplot
locus_scores <- scores(RDA.plate, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores) %>% 
  mutate(type = ifelse(names %in% outliers_plate$Loci, "Outlier", "Non outlier"))

TAB_var <- as.data.frame(scores(RDA.plate, choices=c(1,2), display="bp")) %>% 
  mutate(ID = row.names(scores(RDA.plate, choices=c(1,2), display="bp")) ) #%>% left_join(final.env.names)
TAB_var  # pull the biplot scores

## Biplot of RDA loci and variables scores
gbiplot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  # geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5) +
  ggrepel::geom_label_repel(data = TAB_var, aes(label = ID, x=1.1*RDA1, y=1.1*RDA2), size = 2, max.overlaps = 20
  ) +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.8)),
        strip.background = element_rect(fill = "white"))
gbiplot

ggsave(filename = file.path(filter.ref.path,"01f_Plate", "Plate_biplot_20230419.png"),  # radloci different order than Fis.MeanDepthPerSNP.139317snps.png (graph0.1)
       plot = gbiplot,
       width = 9, height = 7, units = "in", bg = "white")

# Filtration
length(outliers_plate$Loci)

LOC.H06.HW <- read.csv(file.path(filter.ref.path,"01e_HW", "Loc.h06.all.csv"))

LOC.H06.PLATE <- LOC.H06.HW %>% filter(ID %nin% outliers_plate$Loci)

nrow(LOC.H06.PLATE) + length(outliers_plate$Loci) == nrow(LOC.H06.PLATE)

write.csv(LOC.H06.PLATE, file.path(filter.ref.path,"01f_Plate", "Loc.H06.all.PLATE.csv"), 
          row.names = F, quote = F)

vcf.path

cmd6 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"01f_Plate", "Loc.H06.all.PLATE.csv"),
              "--out", file.path(current.wd,filter.ref.path,"01f_Plate", paste0("populations.",LOC.H06.PLATE %>% nrow(),"snps.677ind.H06.all"))  # understand why snps in name are wrongly estimated
)
cmd6

A6 <- system2("vcftools", cmd6, stdout=T, stderr=T)
A6

cat(file = file.path(current.wd,filter.ref.path,"01f_Plate","VCFtools_SnpWhiteList.H06.PLATE.log"),
    "\n", cmd6, "\n",
    A6, # what to put in my file
    append= F, sep = "\n")


if(!file.exists(file.path(filter.ref.path, "01f_Plate", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01f_Plate", ".gitignore")) 
}


## Filtration #7: Sex-linked SNPs ----------------------------------------

list.files(file.path(filter.ref.path, "01f_Plate"))

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path, "01f_Plate", "populations.90117snps.677ind.H06.all.recode.vcf"))

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing sex-linked loci
  select(ID, CHROM, POS) %>% 
  mutate(chromosome = sapply(str_split(CHROM, "_"), `[`,4),
         scaffold = sapply(str_split(CHROM, "_"), `[`,6), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

head(SCAFFOLD.info)
SCAFFOLD.info %>% group_by(chromosome) %>% summarise(N = n())

length(unique(SCAFFOLD.info$RADloc))  # 61159
# RDA on sex

gl.data <- vcfR::vcfR2genlight(vcf.data)

pop.data %>% pull(Sexe_laboratoire) %>% table(useNA = "ifany")

ID.withSex <- pop.data %>% filter(Sexe_laboratoire %in% c("F","M")) %>% pull(ID_GQ)

gl.data.sex <- gl.data[indNames(gl.data) %in% ID.withSex, ]
pop(gl.data.sex) <- data.frame(ID_GQ = indNames(gl.data.sex)) %>% left_join(pop.data) %>% pull(Sexe_laboratoire)

pop(gl.data.sex) %>% table()

library(vegan)

## Full model
RDA.sex <- vegan::rda(tab(gl.data.sex) ~ pop(gl.data.sex))
# tab(gl.data.sex)[1:5,1:5]

anova(RDA.sex)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = tab(gl.data.sex) ~ pop(gl.data.sex))
#           Df Variance      F Pr(>F)    
# Model      1     16.1 1.3092  0.001 ***
# Residual 647   7937.6                  
RsquareAdj(RDA.sex)$adj.r.squared  # 0.0004770025

res.sex <- as.data.frame(scores(RDA.sex, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.sex, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

if(!file.exists(file.path(filter.ref.path, "01g_Sex"))){
  dir.create(file.path(filter.ref.path, "01g_Sex"))
  print(file.path(filter.ref.path, "01g_Sex"))
}

# save(list = c("RDA.sex", "res.sex"),
#      file = file.path(filter.ref.path,"01g_Sex", "SNPs_5x_SexLinked.RData"))
load(file.path(filter.ref.path, "01g_Sex","SNPs_5x_SexLinked.RData"))

res.sex %>% ggplot(aes(x = RDA1, fill = Sexe_laboratoire)) +
  #facet_wrap(~RegionAssesment) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  #geom_point(data = res2, colour = "grey70", cex = 1) +
  geom_histogram() +
  #geom_density()+
  labs(x=c("RDA 1")) +  theme_bw()  

hist(RDA.sex$CCA$v)

sd(RDA.sex$CCA$v)*3

sex.loc <- dimnames(RDA.sex$CCA$v)[[1]][abs(RDA.sex$CCA$v) > 0.01]
sex.loc

div.graph %>% dplyr::mutate(Cat = ifelse(ID %in% outliers_plate$Loci, "Outlier", "Non outlier")) %>% 
  ggplot(aes(x = Hobs, y =  Hexp, col = Cat)) +
  geom_point(alpha = 0.1) +
  facet_grid(~Cat) +
  geom_abline(slope = 1)


# "21171:37:-"   "186296:11:+"  "924530:77:+"  "924531:10:-"  "1395372:9:+"  "1395379:42:-" "1395383:29:-" "1702284:8:+"  "1711667:37:+" "1978242:15:-"
# "1978252:10:+" "2055753:11:+" "2065617:16:+"
# All loci are heterozygous for males (ll 1659-1666)

# Location of sex-linked loci on genome
gl.data.sex$chromosome[gl.data.sex$loc.names %in% sex.loc]

tab(gl.data[, all_of(sex.loc)],NA.method = c("asis")) %>% 
  as.data.frame %>% mutate(ID_GQ = indNames(gl.data)) %>% # head()
  pivot_longer(sex.loc) %>% 
  left_join(pop.data) %>% #head() %>% View()
  left_join(SCAFFOLD.info %>% mutate(CHROM = str_replace(CHROM, "S_20_00703_Chromosome", "Chr") %>% str_replace("scaffold", "sca")) %>% 
              dplyr::select(CHROM, name = ID)) %>% 
  ggplot(aes(x = name, y = ID_GQ, fill = factor(value))) +
  geom_bin2d() +
  facet_grid(Sexe_laboratoire ~ CHROM, space  = "free", scale = "free") + 
  theme_bw() + theme(axis.text.y = element_blank(),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Filtration

HiC.loc <- SCAFFOLD.info %>% dplyr::filter(is.na(chromosome)) %>% pull(ID)
X.loc <- SCAFFOLD.info %>% dplyr::filter(chromosome == "ChromosomeX") %>% pull(ID)

LOC.H06.PLATE <- read.csv(file.path(filter.ref.path,"01f_Plate", "Loc.H06.all.PLATE.csv"))

LOC.H06.SEX <- LOC.H06.PLATE %>% filter(ID %nin% c(HiC.loc, X.loc))

nrow(LOC.H06.SEX) + length(c(HiC.loc, X.loc)) == nrow(LOC.H06.PLATE)

write.csv(LOC.H06.SEX, file.path(filter.ref.path,"01g_Sex", "Loc.H06.all.SEX.csv"), 
          row.names = F, quote = F)

list.files(file.path(here::here(), filter.ref.path, "01f_Plate"))

vcf.path <-  file.path(here::here(), filter.ref.path,"01f_Plate", "populations.90117snps.677ind.H06.all.recode.vcf")

cmd6 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(current.wd,filter.ref.path,"01g_Sex", "Loc.H06.all.SEX.csv"),
              "--out", file.path(current.wd,filter.ref.path,"01g_Sex", paste0("populations.",LOC.H06.SEX %>% nrow(),"snps.677ind.H06.all"))  # understand why snps in name are wrongly estimated
)
cmd6

A6 <- system2("vcftools", cmd6, stdout=T, stderr=T)
A6

cat(file = file.path(current.wd,filter.ref.path,"01g_Sex","VCFtools_SnpWhiteList.H06.SEX.log"),
    "\n", cmd6, "\n",
    A6, # what to put in my file
    append= F, sep = "\n")

if(!file.exists(file.path(filter.ref.path, "01g_Sex", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01g_Sex", ".gitignore")) 
}



### General check: Depth, Position, Homozygosity, and Missing -------------

vcf.path <- file.path(here::here(), filter.ref.path,"01g_Sex", "populations.88433snps.677ind.H06.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing sex-linked loci
  select(ID, CHROM, POS) %>% 
  mutate(chromosome = sapply(str_split(CHROM, "_"), `[`,4),
         scaffold = sapply(str_split(CHROM, "_"), `[`,6), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )
length(unique(SCAFFOLD.info$RADloc))  # 60102

# Verify loci depth after rounds of filtration

gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL)
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))
head(gt.tidy)

vcf.fix <- as.data.frame(vcf.data@fix) %>% mutate(Key = 1:nrow(.))

# group by snps
gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>% 
  left_join(vcf.fix %>% select(Key, ID, CHROM, POS) %>% mutate(POS = as.integer(as.character(POS)))) %>% 
  left_join(lmiss %>% select(CHROM=CHR, POS, F_MISS))

# 'Low' coverage for specific chromosome? It doesn't appear so but for, possibly, some loci around Key 67000
graph6.0 <- gt.key %>% ggplot(aes(x = Key, y = meanDP)) +  # I should see darker spots on a 'column' (Keys) at low meaDP 
  geom_jitter(width = 0, alpha = 0.01) +
  theme_light()

# Estimate Ho/He for remaining loci - first convert to genind format and then estimate Ho and He with summary()
gi.data  <- vcfR::vcfR2genind(vcf.data) 
# Computing diversity (Ho, He)
div <- summary(gi.data)
str(div)

div.graph <- data.frame(ID = names(div$Hobs),
                        Hobs = div$Hobs,
                        Hexp = div$Hexp) %>%
  left_join(gt.key)

div.graph %>% head()

graph6.1 <- div.graph %>% #filter(medianDP < 15) %>% 
  ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.35) +
  scale_colour_distiller(palette = "Spectral") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18, colour = "black"))
graph6.1

ggsave(filename = file.path(filter.ref.path,"01g_Sex", "He_Ho_post_filtering_20230420.png"),
       plot = graph6.1 ,
       width = 10, height = 10, units = "in")

graph6.2 <- div.graph %>% #filter(medianDP < 15) %>% 
  ggplot(aes(x = Hobs, y = Hexp, col = meanDP)) +
  geom_point(alpha = 0.35) +
  scale_color_viridis_c() +
  labs(colour = "Mean Depth") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18, colour = "black"),
        legend.position = c(0.90, 0.125),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13),
        legend.background = element_rect(fill='transparent'))
graph6.2

ggsave(filename = file.path(filter.ref.path,"01g_Sex", "He_Ho_meanDP_post_filtering_20230419.png"),
       plot = graph6.2 ,
       width = 10, height = 10, units = "in")

graph6.3 <- div.graph %>% #filter(medianDP < 15) %>% 
  ggplot(aes(x = Hobs, y = Hexp, col = medianDP)) +
  geom_point(alpha = 0.35) +
  scale_color_viridis_c() +
  labs(colour = "Median Depth") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18, colour = "black"),
        legend.position = c(0.90, 0.125),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13),
        legend.background = element_rect(fill='transparent'))
graph6.3

ggsave(filename = file.path(filter.ref.path,"01g_Sex", "He_Ho_medianDP_post_filtering_20230420.png"),
       plot = graph6.3 ,
       width = 10, height = 10, units = "in")

graph6.4 <- div.graph %>% filter(meanDP > 19) %>%  # (mean(div.graph$meanDP) = 18.92) from 74k loci to 31.9k loci
  ggplot(aes(x = Hobs, y = Hexp, col = meanDP)) +
  geom_point(alpha = 0.35) +
  scale_color_viridis_c() +
  labs(colour = "Median Depth") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18, colour = "black"),
        legend.position = c(0.90, 0.125),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13),
        legend.background = element_rect(fill='transparent'))
graph6.4

ggsave(filename = file.path(filter.ref.path,"01g_Sex", "He_Ho_byMeanDP19_post_filtering_20230419.png"),
       plot = graph6.4 ,
       width = 10, height = 10, units = "in")

graph6.5 <- div.graph %>% #filter(F_MISS > 15) %>% 
  ggplot(aes(x = Hobs, y = Hexp, col = F_MISS)) +
  geom_point(alpha = 0.35) +
  scale_color_viridis_c() +
  labs(colour = "Proportion missing alleles") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 18, colour = "black"),
        legend.position = c(0.85, 0.125),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13),
        legend.background = element_rect(fill='transparent'))
graph6.5

ggsave(filename = file.path(filter.ref.path,"01g_Sex", "He_Ho_byFMISS_post_filtering_230203.png"),
       plot = graph6.5 ,
       width = 10, height = 10, units = "in")

# div.graph %>% filter(Key > 50000) %>%   # filter by position in alignment
#   ggplot(aes(x = Hobs, y = Hexp, col = F_MISS)) +
#   geom_point(alpha = 0.35) +
#   scale_color_viridis_c() +
#   #labs(colour = "Proportion missing alleles") +
#   theme_bw(base_size = 20, base_family = "Helvetica") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 18, colour = "black"),
#         #legend.position = c(0.85, 0.125),
#         #legend.title = element_text(size = 14), 
#         #legend.text = element_text(size = 13),
#         #legend.background = element_rect(fill='transparent')
#   )



## Filtration #8: Relatedness ---------------------------------------------

current.wd <- here::here()
file.path(here::here(),filter.ref.path,"01g_Sex") %>% list.files()

# In VCF tool, check for duplicated ind.

cmd7a <- paste("--vcf", file.path(here::here(), filter.ref.path,"01g_Sex", "populations.88433snps.677ind.H06.all.recode.vcf"), 
               #sub_indv,
               # --depth,
               "--relatedness2"
)

if(!file.exists(file.path(filter.ref.path, "01h_Relatedness"))){
  dir.create(file.path(filter.ref.path, "01h_Relatedness"))
  print(file.path(filter.ref.path, "01h_Relatedness"))
}

setwd(file.path(filter.ref.path,"01h_Relatedness")) 

A7a <- system2("vcftools", cmd7a, stdout=T, stderr=T)
A7a

cat(file = file.path("VCFtools.Relatedness.log"),
    "\n", cmd7a, "\n",
    A7a, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

related2.all <- read.delim(file.path(filter.ref.path,"01h_Relatedness","out.relatedness2"), header = T, sep = "\t" )
related2.all %>% head()

length(related2.all$INDV1 %>% unique)

# Remove identical individuals

related2.ok <- related2.all %>% mutate(iden = ifelse(INDV1 == INDV2, 1, 0)) %>%  # identify same ID in INDIV1 and INDIV2
  filter(iden == 0) %>%  # Remove identical IDs
  left_join(pop.data, by = c("INDV1" = "ID_GQ"))  # attach metadata

related2.ok %>% filter(RELATEDNESS_PHI >=1/4) %>% arrange(desc(RELATEDNESS_PHI)) %>% View()

ID.relatedness.checklist <-  related2.ok %>% filter(RELATEDNESS_PHI >=1/4) %>% arrange(desc(RELATEDNESS_PHI)) %>% 
  select(INDV1, INDV2, RELATEDNESS_PHI, No_plaque_envoi, No_puits_envoi, Numero_unique_specimen, Sexe_laboratoire,
         Region_echantillonnage, Lieu_echantillonnage, Latitude_echantillonnage_DD, Longitude_echantillonnage_DD,
         Annee_echantillonnage, Mois_echantillonnage, Jour_echantillonnage)

write_csv(ID.relatedness.checklist,
          file.path(filter.ref.path,"01h_Relatedness","ID.relatedness.checklist.csv"))

rows <- c("09", "10")
ID.clones <- ID.relatedness.checklist[!(ID.relatedness.checklist$No_plaque_envoi %in% "PE_20_00003" & 
                                        grepl(paste(rows, collapse = '|'), ID.relatedness.checklist$No_puits_envoi, ignore.case = T)),]


# Verify is high relatedness is caused by high proportion of NAs

cmd7b <- paste("--vcf", file.path(current.wd, filter.ref.path, "01g_Sex", "populations.88433snps.677ind.H06.all.recode.vcf"), 
               #"--out",  "MAX2.NArm",
               #" --max-alleles", 2,
               #"--max-missing", 0.80,
               #"--missing-site",
               "--missing-indv"
               #"--kept-sites"
               #"--recode"
)
cmd7b

setwd(file.path(filter.ref.path, "01h_Relatedness")) 

A7b <- system2("vcftools", cmd7b, stdout=T, stderr=T)
A7b

cat(file = "populations.filt_GLOBAL_5x_Individuals_wMissing.log",
    "\n", cmd7b, "\n",
    A7b, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

imiss <- read.delim(file.path(filter.ref.path, "01h_Relatedness", "out.imiss"), skip=0, sep = "\t", header = T )
imiss %>% head()
imiss %>% filter(F_MISS >.10)  # of interest: S_201060, S_20_0713, S_20_0731, S_20_1047

imiss %>% left_join(pop.data %>% select(INDV = ID_GQ, Espece)) %>%  ggplot(aes(x=F_MISS, fill = Espece)) + geom_histogram()

imiss %>% left_join(pop.data, by = c("INDV" = "ID_GQ")) %>%
  ggplot(aes(x = Region_echantillonnage, y = F_MISS)) + 
  geom_boxplot(col = "black", fill = NA) +
  geom_jitter(height = 0, alpha = 0.2, col = "blue") +
  geom_hline(yintercept = 0.30, lty = "dashed", col = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

# List of ID to remove

ID.to.remove <- c(ID.relatedness.checklist$INDV1, ID.relatedness.checklist$INDV2) %>% unique()  # this lists (indiscriminately) all IDs with phi > 0.25 (parent-offspring)
# In MOBELS, some extracts belong to the same ID and the high phi is not caused by cross-contamination - remove only one of the two
# When the high phi is caused by cross contamination, remove both IDs (e.g. PE_20_00003, rows 9 and 10)
ID.to.remove2 <- c("S_20_01631","S_20_02444","S_20_03119","S_20_03575","S_20_00796","S_20_03524","S_20_00762","S_20_02381","S_20_03094","S_20_02913",
                   "S_20_02716","S_20_03616","S_20_00731","S_20_00759","S_20_01052","S_20_01053","S_20_00743","S_20_01102","S_20_01112","S_20_01006",
                   "S_20_03038","S_20_01059","S_20_01060")

related2.ok %>% filter(RELATEDNESS_PHI >=1/4) %>% arrange(desc(RELATEDNESS_PHI)) %>% 
  left_join(imiss, by = c("INDV1" = "INDV")) %>% select(c(1,2,7,10:11,30)) %>% View(.)

# ID.to.remove <- c(ID.relatedness.checklist$INDV1, ID.relatedness.checklist$INDV2) %>% unique()  # unsure why am I creating this since I'm not using it anymore as it contains ID that should remain in dataset
ID.plate3 <- pop.data %>% filter(No_plaque_envoi == "PE_20_00003" , No_puits_envoi %in% c(paste0(LETTERS[1:8], "09"), paste0(LETTERS[1:8], "10"))) %>% 
  pull(ID_GQ)  # list of individuals in plate 3 from rows 9 and 10

ID.relatedness.checklist %>% filter(No_plaque_envoi == "PE_20_00003") %>% pull(No_puits_envoi) %>% unique() %>% sort()

related2.ok %>% arrange(desc(RELATEDNESS_PHI)) %>% head()

c(related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% pull(INDV1),
  related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% pull(INDV2)
) %>% table()

graph.relatedness <- related2.ok %>% 
  #filter(INDV1 %nin% c(ID.to.remove,ID.plate6)  ,
  #       INDV2 %nin% c(ID.to.remove,ID.plate6) ) %>% 
  ggplot(aes(x =  RELATEDNESS_PHI)) +
  labs(x = "Relatedness coefficient", y = "N observations", title = "Relatedness coefficient distribution")+
  geom_histogram(bins = 100) +
  geom_vline(xintercept = c(1/2, 1/4, 1/8, 1/16), lty = "dashed")+
  #facet_wrap(~No_plaque) +
  annotate("text", 
           x = c(1/2-0.005, 1/4-0.005, 1/8-0.005, 1/16-0.005), y = 1000, 
           label = c("Individual-self", "Siblings / Parent-offspring", "Half-siblings / Grandparent-grandchild", "First cousins"), 
           angle = 90, hjust = 0, vjust = 0) +
  theme_bw()  + #
  ggforce::facet_zoom(xlim = c(0.125, 0.5), ylim = c(0, 400), zoom.size = 1, horizontal = FALSE)
#  theme(axis.text.x = element_text(angle = 60, hjust = 1))

graph.relatedness 

ggsave(filename = file.path(filter.ref.path,"01h_Relatedness", "Relatedness230203.png"),
       plot = graph.relatedness ,
       width = 10, height = 10, units = "in")

# Remove duplicates specimens (same ID, cross-contamination, etc)

good.ID <- related2.ok %>% 
  filter(INDV1 %nin% unique(c(ID.to.remove2,ID.plate3)),
         INDV2 %nin% unique(c(ID.to.remove2,ID.plate3))) %>% pull(INDV1, INDV2) %>% unique() 

# good.ID  %>% length()  # 638 IDs

# pop.data %>% dplyr::filter(ID_GQ %in% good.ID) %>% dplyr::select(ID_GQ, Numero_unique_specimen, Region_echantillonnage) %>% write_csv("Final.list.Luca.20230420.csv")
# 
# 
# TREV <- c("S_20_01788", 	"S_20_03480", 	"S_20_01766", 	".S_20_03664", 	"S_20_01763", 	"S_20_00622", 	"S_20_00710", 	"S_20_00711", 	"S_20_03448", 	"S_20_03453", 	"S_20_00661", 	"S_20_00691", 	"S_20_00694", 	"S_20_00776", 	"S_20_00785", 	"S_20_00702", 	"S_20_00621", 	"S_20_00693", 	"S_20_03515", 	"S_20_01067", 	"S_20_01019", 	"S_20_00751")
# TREV
# 
# TREV[TREV %in% good.ID]

list.files(file.path(here::here(), filter.ref.path,"01g_Sex"))

vcf.path <-  file.path(here::here(), filter.ref.path,"01g_Sex", "populations.88433snps.677ind.H06.all.recode.vcf")

cmd7c <- paste("--vcf", vcf.path, 
               "--recode",
               paste("--indv", good.ID, collapse = " "),
               "--out",  vcf.path %>% str_replace("01g_Sex", "01h_Relatedness") %>% 
                 str_replace("677ind", paste0(length(good.ID), "ind")) %>% 
                 str_remove(".recode.vcf")
)
cmd7c

A7c <- system2("vcftools", cmd7c, stdout=T, stderr=T)
A7c %>% tail()

cat(file = file.path(filter.ref.path, "01h_Relatedness","VCFtools_RemoveRelatedness.log"),
    "\n", cmd7c, "\n",
    A7c, # what to put in my file
    append= F, sep = "\n")


# Add gitignore if necessary

if(!file.exists(file.path(filter.ref.path, "01h_Relatedness", ".gitignore")) ){
  cat("out.relatedness2", "*.vcf", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01h_Relatedness", ".gitignore")) 
}



# Final dataset -----------------------------------------------------------

list.files(file.path(filter.ref.path, "01h_Relatedness"))

vcf.path <- file.path(filter.ref.path, "/01h_Relatedness/populations.88433snps.638ind.H06.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

#gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# Transform it into a data frame
SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing sex-linked loci
  select(ID, CHROM, POS) %>% 
  mutate(chromosome = sapply(str_split(CHROM, "_"), `[`,4),
         scaffold = sapply(str_split(CHROM, "_"), `[`,6), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

# SCAFFOLD.info %>% group_by(CHROM) %>% summarise(N = n()) %>% print (n = 500)


#ID.final <-  data.frame(ID_GQ = indNames(gi.data)) %>% 
#  left_join(pop.data) %>% 
#  filter(Gen_ZONE_FG %nin% c("SFA-6-3", "SFA-8-4", "SFA-12-2", "SFA-12-3", "SFA-15-1")) %>% 
#  pull(ID_GQ)

#length(ID.final)

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

# Check that the order was kept
table(locNames(gi.data) == SCAFFOLD.info$ID)

# ADD NA and MAF info to select the "BEST" SNP / RAD loci

SCAFFOLD.info$NNA <- na.gi.count(gi.data)
SCAFFOLD.info$MAF <- adegenet::minorAllele(gi.data)

SCAFFOLD.info %>% group_by(RADloc) %>% summarise(N = n()) %>% arrange(desc(N)) %>% head()

SCAFFOLD.info.final <- SCAFFOLD.info %>% mutate(RADloc = as.numeric(RADloc)) %>% 
  arrange(RADloc, desc(round(MAF, 1)), round(NNA,2)) %>% 
  distinct(RADloc, .keep_all = T)

SCAFFOLD.info$MAF %>% hist()
SCAFFOLD.info.final$MAF %>% hist()

# # Function to create a list of loci, from a genind object
# 
# filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
#   # Create vectors for each loci
#   MAF.res <- adegenet::minorAllele(gi)
#   NA.res  <- na.gi.count(gi)
#   
#   # Filter by threshold
#   MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
#   cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
#   
#   NA.loc <- names(NA.res[NA.res <= NA.trs])
#   cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
#   
#   # LOCI with both conditions
#   LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
#   LOCI.res %>% length()
#   
#   cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
#   
#   return(LOCI.res)
# }
# 

if(!file.exists(file.path(filter.ref.path, "01i_UniqueFinal"))){
  dir.create(file.path(filter.ref.path, "01i_UniqueFinal"))
  print(file.path(filter.ref.path, "01i_UniqueFinal"))
}

write.csv(SCAFFOLD.info.final %>% select(ID), file.path(filter.ref.path,"01i_UniqueFinal", "Loc.H06.DP.unique.final.csv"), 
          row.names = F, quote = F)

write.csv(SCAFFOLD.info.final, file.path(filter.ref.path,"01i_UniqueFinal", "Scaffold.info.H06.DP.unique.final.csv"), 
          row.names = F, quote = F)


# CREATE VCF WITH UNIQUE RAD LOCI

vcf.path

cmd10 <- paste("--vcf", vcf.path, 
               "--recode",
               "--snps", file.path(here::here(),filter.ref.path, "01i_UniqueFinal", "Loc.H06.DP.unique.final.csv"),
               
               "--out", file.path(here::here(),filter.ref.path, "01i_UniqueFinal", paste0("populations.", SCAFFOLD.info.final %>% nrow(),"snps.","638ind.single.final"))
)
cmd10

A10 <- system2("vcftools", cmd10, stdout=T, stderr=T)
tail(A10)

cat(file =  file.path(filter.ref.path, "01i_UniqueFinal", "VCF.filter.log"),
    "\n", cmd10, "\n",
    A10, # what to put in my file
    append= F, sep = "\n")


# Reload, and save as Rdata
list.files(file.path(here::here(),filter.ref.path, "01i_UniqueFinal"))

vcf.path <- file.path(here::here(),filter.ref.path, "01i_UniqueFinal/populations.60102snps.638ind.single.final.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing sex-linked loci
  select(ID, CHROM, POS) %>% 
  mutate(chromosome = sapply(str_split(CHROM, "_"), `[`,4),
         scaffold = sapply(str_split(CHROM, "_"), `[`,6), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

length(unique(SCAFFOLD.info$RADloc))  # 60102

gl.final <- vcfR::vcfR2genlight(vcf.data) 
gi.final <- vcfR::vcfR2genind(vcf.data) 

gl.final 
gi.final

pop(gl.final) <- data.frame(ID_GQ = indNames(gl.final)) %>% 
  left_join(pop.data) %>% 
  pull(Region_echantillonnage)

table(pop(gl.final))

pop(gi.final) <- data.frame(ID_GQ = indNames(gi.final)) %>% 
  left_join(pop.data) %>% 
  pull(Region_echantillonnage)

save(list = c("gl.final", "gi.final"),
     file = file.path(here::here(),filter.ref.path, "01i_UniqueFinal/populations.60102snps.638ind.H06.DP.single.final.recode.vcf.adegenet.Rdata"))

if(!file.exists(file.path(filter.ref.path, "01i_UniqueFinal", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "01i_UniqueFinal", ".gitignore")) 
}
