# Info --------------------------------------------------------------------
#
# Overview: Testing ADMIXTURE on the ddRAD 27,000 SNPs dataset
# 
# Author: Audrey Bourret
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-04-05
#


# Library -----------------------------------------------------------------

library(tidyverse)


# Dataset -----------------------------------------------------------------

pop.data <- read.csv(file = "./00_Data/02_Dataset/Beluga_ddRAD.csv")


# Admixture ---------------------------------------------------------------

# Original SNP

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD/populations.27352snps.639ind.final.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)


for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD/02b_Admixture/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               "-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("test.k",k, ".log"),
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
  temp <- readLines(file.path("./02_Results/01_ddRAD/02b_Admixture/", paste0("test.k",k, ".log")))
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
gg.CV

ggsave(filename = file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", "Admixture.test.CV.png"), 
       plot = gg.CV,
       height = 3.5, width = 4, units = "in")   


k <- 6

Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", paste0("populations.27352snps.639ind.final.recode.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", paste0("populations.27352snps.639ind.final.recode.",3,".Q")))
Q.k4.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", paste0("populations.27352snps.639ind.final.recode.",4,".Q")))
Q.k5.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", paste0("populations.27352snps.639ind.final.recode.",5,".Q")))
Q.k6.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/", paste0("populations.27352snps.639ind.final.recode.",6,".Q")))


Q.res <- bind_rows(cbind(fam$V1, Q.k6.res, K = 6),
                   cbind(fam$V1, Q.k5.res, K = 5),
                   cbind(fam$V1, Q.k4.res, K = 4),
                   cbind(fam$V1, Q.k3.res, K = 3),
                   cbind(fam$V1, Q.k2.res, K = 2))


head(Q.res)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res) <- c("ID_GQ", paste0("Q", 1:6), "K")

#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str.all <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>% 
#  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
#  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  left_join(pop.data) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Longitude_echantillonnage_DD, FUN = function(x) max(as.numeric(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
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


ggsave(filename = "02_Results/01_ddRAD/02b_Admixture/Admixture_k2_to_k6.png", plot = gg.str.all,
       width = 12, height = 6, unit = "in")

gg.str.k4 <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
    dplyr::filter(K == 6) %>% 
  left_join(pop.data) %>% 
  mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Season, FUN = function(x) max(levels(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(. ~ Region2 + Season, space = "free", scale = "free") +
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
gg.str.k4

ggsave(filename = "02_Results/01_ddRAD/02b_Admixture/Admixture_k4_season.png", plot = gg.str.k4,
       width = 12, height = 3, unit = "in")



# Trying maps -------------------------------------------------------------

library(terra)
library(ggspatial)
library(tidyterra)

admin <- terra::vect(rnaturalearth::ne_countries(scale = "medium", continent = "north america", returnclass	  = "sf"))


pop.df <- pop.data %>% 
  mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
Region = ifelse(Region == "LON", "JAM", Region)) %>% 
  group_by(Region2, Region) %>% summarise(Lon = mean(Lon,  na.rm = TRUE), Lat = mean(Lat,  na.rm = TRUE),
                                          N = n()) 

pop.df[which(pop.df$Region == "REB"), c("Lon")] <- c(-86.3) #, 66.6)
pop.df[which(pop.df$Region == "REB"), c("Lat")] <- c(66.6)
pop.df[which(pop.df$Region == "SWH"), c("Lon")] <- c(-92.24) #, 66.6)
pop.df[which(pop.df$Region == "SWH"), c("Lat")] <- c(57.2)

pop.df


pop <- terra::vect(pop.df,
                   geom = c("Lon", "Lat"),  crs = "+proj=longlat")

plot(crop(admin, ext(pop)))
plot(pop, add = T)

gg.env <- ggplot(pop) + 
  geom_sf(data = admin, fill="gray95", size=0.5, alpha = 1) +
  
  geom_sf(data = pop, alpha = .8, aes(col = Region2 , size = N )) +  
  #scale_color_viridis_c(name = "Depth (m)", direction = -1, option = "plasma")+
  # Map limits
  #scale_shape_manual(values = c(16,15,17), name = "Year") +
  coord_sf(xlim = c(-2000000, 1000000), ylim = c(500000, 3000000), crs = sf::st_crs("EPSG:6622")) +
  # Others
  #facet_grid(. ~ pred.pop) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Beluga ddRAD") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom")  
gg.env





Q.pie <-  Q.res %>% pivot_longer(cols =  paste0("Q", 1:6), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
  dplyr::filter(K == 6) %>% 
  left_join(pop.data) %>% 
  mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
         Region = ifelse(Region == "LON", "JAM", Region),
         Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter"))))) %>% 
  select(-c(Lat, Lon)) %>% 
  left_join(pop.df) %>% #View()
  group_by( Region2, Region,  Lat, Lon, Season, Group) %>% 
  summarise(value = mean(Q)# ,
            #Lat = mean(Lat),
            #Long = mean(Lon)
            )           %>% 
left_join(pop.data %>%   mutate(Region2, factor(Region2, levels = c("SLE","CSB","FRB","UNG","SHS","NHS","NEH","EHB","BEL","JAM","WHB")),
                                Region = ifelse(Region == "LON", "JAM", Region),
                                Season = ifelse(is.na(Month), "Unknown",
                                                ifelse(Month %in% c(7,8), "Summer",
                                                       ifelse(Month > 3 & Month < 7, "Spring",
                                                              ifelse(Month > 8 & Month < 12, "Fall",
                                                                     "Winter"))))) %>%  group_by(Region, Season) %>% summarise(N = n())) %>% 
  dplyr::filter(Season != "Unknown")
  
  #  mutate(N = sum(value, na.rm = T))
Q.pie %>% View()

#GSL.ID %>%  ggplot(aes(x =as.numeric(Longueur_mm) , fill =  factor(Annee_echantillonnage))) + geom_histogram() + facet_grid(pred.pop ~ .)


library(scatterpie)
gg.map <- ggplot() + 
  geom_sf(data = admin, fill=NA, cex=0.1) +
  geom_scatterpie(aes(x=Lon, y=Lat, group = Region, r = 1.5), 
                  data = Q.pie, cols = "Group",long_format = T) +
  coord_sf(xlim = c(-100, -60), ylim = c(45, 70)) +
  scale_fill_brewer(palette = "Set1") + 
  facet_grid(~Season) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

gg.map

ggpubr::ggarrange( gg.str.k4 + theme(legend.position = "none"), gg.map, nrow = 2, common.legend = F
                 

                   )

# Admixture without SLE ---------------------------------------------------
file.path(here::here(), "./00_Data/03_ddRAD/populations.27352snps.639ind.final.recode.bed")


# MEN


ID.noSLE <- Q.res %>% dplyr::filter(K == 2) %>% left_join(pop.data) %>% dplyr::filter(Region2 != "SLE") %>% pull(ID_GQ)
length(ID.noSLE)

cmd <- paste("--vcf", file.path(file.path("./00_Data/03_ddRAD/populations.27352snps.639ind.final.recode.vcf")), 
             "--recode",
             paste("--indv",ID.noSLE, collapse = " "),
             "--out", file.path("./00_Data/03_ddRAD/01_Subset_dataset/populations_noSLE.27352snps.616ind.final")
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)


# Convert as plink tped too for pcadapt

vcf.path <- file.path("./00_Data/03_ddRAD/01_Subset_dataset/populations_noSLE.27352snps.616ind.final.recode.vcf")
file.exists(vcf.path)

cmd1 <- paste("--vcf", vcf.path, 
              #"--recode",
              "--plink-tped",
              "--out",  vcf.path %>% str_remove(".vcf"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", "./00_Data/03_ddRAD/01_Subset_dataset/populations_noSLE.27352snps.616ind.final.recode.tfam", 
               "--tped", "./00_Data/03_ddRAD/01_Subset_dataset/populations_noSLE.27352snps.616ind.final.recode.tped", 
               "--make-bed", 
               "--out", vcf.path %>% str_remove(".vcf")
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a


# Original SNP

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD/01_Subset_dataset/populations_noSLE.27352snps.616ind.final.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)


for(k in 1:10){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD/02b_Admixture/01_NoSLE") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("test.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.res <- data.frame(k = 1:6,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("test.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.res[i, "CV"] <- CV
  
}

CV.res$CV <- as.numeric(as.character(CV.res$CV))

CV.res %>% arrange(CV)

plot(CV.res$CV)

gg.CV <- CV.res %>% mutate(color = ifelse(k == 3, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
gg.CV

ggsave(filename = file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE/", "Admixture.test.CV.png"), 
       plot = gg.CV,
       height = 3.5, width = 4, units = "in")   


k <- 5

Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",3,".Q")))
Q.k4.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",4,".Q")))
Q.k5.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",5,".Q")))


Q.res <- bind_rows(#cbind(fam$V1, Q.k6.res, K = 6),
                   cbind(fam$V1, Q.k5.res, K = 5),
                   cbind(fam$V1, Q.k4.res, K = 4),
                   cbind(fam$V1, Q.k3.res, K = 3),
                   cbind(fam$V1, Q.k2.res, K = 2))


head(Q.res)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res) <- c("ID_GQ", paste0("Q", 1:5), "K")

#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str.all <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:5), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
    dplyr::filter(ID_GQ %in% ID.largeGroup) %>% 
  left_join(pop.data) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Longitude_echantillonnage_DD, FUN = function(x) max(as.numeric(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
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


ggsave(filename = "02_Results/01_ddRAD/02b_Admixture/01_NoSLE/Admixture_k2_to_k5.png", plot = gg.str.all,
       width = 12, height = 6, unit = "in")





# Admixture large undefined group ---------------------------------------------------

# MEN


ID.largeGroup <- Q.res %>% dplyr::filter(K == 4, Q2 >= 0.50) %>% pull(ID_GQ)
length(ID.largeGroup)

cmd <- paste("--vcf", file.path(file.path("./00_Data/03_ddRAD/populations.27352snps.639ind.final.recode.vcf")), 
             "--recode",
             paste("--indv",ID.largeGroup, collapse = " "),
             "--out", file.path("./00_Data/03_ddRAD/01_Subset_dataset/populations_largeGroup.27352snps.494ind.final")
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)

# Convert as plink tped too for pcadapt

vcf.path <- file.path("./00_Data/03_ddRAD/01_Subset_dataset/populations_largeGroup.27352snps.494ind.final.recode.vcf")
file.exists(vcf.path)

cmd1 <- paste("--vcf", vcf.path, 
              #"--recode",
              "--plink-tped",
              "--out",  vcf.path %>% str_remove(".vcf"))


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)

cmd2a <- paste("--tfam", "./00_Data/03_ddRAD/01_Subset_dataset/populations_largeGroup.27352snps.494ind.final.recode.tfam", 
               "--tped", "./00_Data/03_ddRAD/01_Subset_dataset/populations_largeGroup.27352snps.494ind.final.recode.tped", 
               "--make-bed", 
               "--out", vcf.path %>% str_remove(".vcf")
               
)

A2a <- system2("plink", cmd2a, stdout=T, stderr=T)
A2a


# Original SNP

bed.file <- file.path(here::here(), "./00_Data/03_ddRAD/01_Subset_dataset/populations_largeGroup.27352snps.494ind.final.recode.bed")
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)


for(k in 1:3){
  
  print(k)  
  
  setwd(file.path(here::here(), "/02_Results/01_ddRAD/02b_Admixture/02_LargeGroup") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               #"-B999",
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("n494.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.res <- data.frame(k = 1:3,
                     CV = NA,
                     stringsAsFactors = F)


for(i in 1:nrow(CV.res)){
  # Which k
  k <- CV.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_ddRAD/02b_Admixture/02_LargeGroup", paste0("n494.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.res[i, "CV"] <- CV
  
}

CV.res$CV <- as.numeric(as.character(CV.res$CV))

CV.res %>% arrange(CV)

plot(CV.res$CV)

gg.CV <- CV.res %>% mutate(color = ifelse(k == 3, "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
gg.CV

ggsave(filename = file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE/", "Admixture.test.CV.png"), 
       plot = gg.CV,
       height = 3.5, width = 4, units = "in")   


k <- 3

Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/02_LargeGroup", paste0("populations_largeGroup.27352snps.358ind.final.recode.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/02_LargeGroup", paste0("populations_largeGroup.27352snps.358ind.final.recode.",3,".Q")))
#Q.k4.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",4,".Q")))
#Q.k5.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/01_NoSLE", paste0("populations_noSLE.27352snps.616ind.final.recode.",5,".Q")))

Q.k2.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/02_LargeGroup", paste0("populations_largeGroup.27352snps.494ind.final.recode.",2,".Q")))
Q.k3.res <-  read.table(file.path(here::here(), "02_Results/01_ddRAD/02b_Admixture/02_LargeGroup", paste0("populations_largeGroup.27352snps.494ind.final.recode.",3,".Q")))



Q.res <- bind_rows(#cbind(fam$V1, Q.k6.res, K = 6),
  #cbind(fam$V1, Q.k5.res, K = 5),
  #cbind(fam$V1, Q.k4.res, K = 4),
  cbind(fam$V1, Q.k3.res, K = 3),
  cbind(fam$V1, Q.k2.res, K = 2))


head(Q.res)

#Q.fas.res <- cbind(fam.fas$V1, Q.fas.res)

names(Q.res) <- c("ID_GQ", paste0("Q", 1:3), "K")

#reorder(ID, Qvalue, FUN = function(x) min(x))

gg.str.all <- Q.res %>% pivot_longer(cols =  paste0("Q", 1:3), names_to = "Group", values_to = "Q") %>% 
  #  mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
  #  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  left_join(pop.data) %>% 
  #ggplot(aes(x = reorder(ID_GQ, Longitude_echantillonnage_DD, FUN = function(x) max(as.numeric(x))), y = Q, fill = Group)) + 
  ggplot(aes(x = ID_GQ, y = Q, fill = Group)) + 
  
  geom_col() +
  #facet_grid(. ~Lieu_echantillonnage + Mois_echantillonnage, space = "free", scale = "free") +
  facet_grid(K ~ Region , space = "free", scale = "free") +
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


ggsave(filename = "02_Results/01_ddRAD/02b_Admixture/02_LargeGroup/Admixture_k2_to_k3.png", plot = gg.str.all,
       width = 12, height = 6, unit = "in")
