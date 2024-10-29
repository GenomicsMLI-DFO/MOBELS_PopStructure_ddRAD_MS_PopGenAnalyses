# Info --------------------------------------------------------------------
#
# Overview: Merged results (PCA, Admixture, Fst) from both ddRAD and lcWGS datasets
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory
# Location: Maurice Lamontagne Institute
# Date: 2023-07-20
#
# Note: Panel B of Figure 4 is made here
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

# Functions
"%nin%" <- Negate("%in%")




# Data --------------------------------------------------------------------

## ddRAD ------------------------------------------------------------------

# Metadata
rad.pop.data <- read.csv(file = "./00_Data/02_Dataset/Beluga_ddRAD.csv") %>% select(1,3:4,9:25) %>% 
  mutate(ID_GQ = gsub("_rep", "", ID_GQ))

# PCA
rad.PCA <- read.csv("./02_Results/01_ddRAD_Bringloe/01_PCA/PCA.neutral.MAF05NA05_n638.csv", stringsAsFactors = F) %>% select(1:9) %>% 
  mutate(ID = gsub("_rep", "", ID))

# Admixture
rad.admx <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_neutral_K5.csv", stringsAsFactors = F) %>% select(1,11:16) %>% 
  mutate(ID_GQ = gsub("_rep", "", ID_GQ))
# rad.admx.pan <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_PAN_K3.csv", stringsAsFactors = F) %>% select(1:4,12) %>% 
#   mutate(ID_GQ = gsub("_rep", "", ID_GQ))
# rad.admx.jam <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_JAM_K2.csv", stringsAsFactors = F) %>% select(1:3,12) %>% 
#   mutate(ID_GQ = gsub("_rep", "", ID_GQ))

# Fst
rad.fst <- read.csv("./02_Results/01_ddRAD_Bringloe/03_FST/FST_ddRAD_neutral.csv", stringsAsFactors = F)
# rad.fst.pan <- read.csv("./02_Results/01_ddRAD_Bringloe/03_FST/", stringsAsFactors = F)
# rad.fst.jam <- read.csv("./02_Results/01_ddRAD_Bringloe/03_FST/", stringsAsFactors = F)

# Merge ddRAD results together
d.rad <- rad.pop.data %>% filter(ID_GQ %in% rad.PCA$ID) %>% filter(!duplicated(ID_GQ)) %>%  # remove unused ID and duplicates
  left_join(rad.PCA, by = c("ID_GQ"="ID")) %>% 
  left_join(rad.admx, by = "ID_GQ") #%>% 
  # left_join(rad.admx.pan, by = "ID_GQ") %>% 
  # left_join(rad.admx.jam, by = "ID_GQ")
colnames(d.rad)[c(2,8:9,18,21:28,34)] <- c("ID_ext","Sex.visual","Sex.qPCR","Haplotype","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","Membership") 
# d.rad <- d.rad %>% mutate(Membership = ifelse(is.na(Membership.jam) & is.na(Membership.pan), Membership.all,
#                                               ifelse(is.na(Membership.jam) & !is.na(Membership.pan), Membership.pan,
#                                                      ifelse(!is.na(Membership.jam), Membership.jam,
#                                                             NA))), 
#                           Marker = "ddRAD")
d.rad %>% group_by(Membership) %>% summarise(N = n())
d.rad <- d.rad %>% mutate(Membership = gsub("CSB","CS",Membership),
                          Marker = "ddRAD")


## lcWGS ------------------------------------------------------------------

# Metadata
wgs.pop.data <- read.csv(file = "./00_Data/02_Dataset/MOBELS_WGS_230717.csv", stringsAsFactors = F) %>% filter(!duplicated(Numero_unique_specimen)) %>% select(1,28,9:12,14:15,25,
                                                                                                                                                    16:23,29)
colnames(wgs.pop.data) <- c("ID","ID_ext","ID_rcpt","Common_name","Genus","Species","Age","Sex.visual","Sex.qPCR","Region","Location","Lat","Lon",
                            "Year","Month","Day","Community","Haplotype")
wgs.pop.data <- wgs.pop.data %>% mutate(Region1 = ifelse(Region %in% "CSB", "CSB",
                                                  ifelse(Region %in% "FRB", "FRB",
                                                  ifelse(Region %in% c("JAM","LON"), "JAM",
                                                  ifelse(Region %in% "NEH", "NEH",
                                                  ifelse(Region %in% "NHB", "NHB",
                                                  ifelse(Region %in% "NHS", "NHS",
                                                  ifelse(Region %in% "RES", "RES",
                                                  ifelse(Region %in% "SAN", "BEL",
                                                  ifelse(Region %in% "SEH", "SEH",
                                                  ifelse(Region %in% "SHS", "SHS",
                                                  ifelse(Region %in% "SLE", "SLE",
                                                  ifelse(Region %in% "UNG", "UNG",
                                                  ifelse(Region %in% "NWH", "NWH",
                                                  ifelse(Region %in% "SWH", "SWH",
                                                  NA))))))))))))))) %>% 
  mutate(Region2 = ifelse(Region %in% "CSB", "CSB",
                   ifelse(Region %in% "FRB", "FRB",
                   ifelse(Region %in% c("JAM","LON"), "JAM",
                   ifelse(Region %in% "NEH", "NEH",
                   ifelse(Region %in% c("NHB","NWH","SWH"), "WHB",
                   ifelse(Region %in% "NHS", "NHS",
                   ifelse(Region %in% "RES", "RES",
                   ifelse(Region %in% "SAN", "BEL",
                   ifelse(Region %in% "SEH", "SEH",
                   ifelse(Region %in% "SHS", "SHS",
                   ifelse(Region %in% "SLE", "SLE",
                   ifelse(Region %in% "UNG", "UNG",
                   NA)))))))))))))

# PCA
wgs.PCA <- read.csv("./02_Results/01_ddRAD_Bringloe/01_PCA/MOBELS-S_20_00703-Beluga-lcWGS-PCA-17vii23.csv", stringsAsFactors = F) %>% select(1,3:10)

# Admixture
wgs.admx <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/MOBELS-S_20_00703-Beluga-lcWGS-Metadata+K6+Het-9v23.csv", stringsAsFactors = F) %>% select(2,14:20)
# wgs.admx.pan <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/MOBELS-S_20_00703-Beluga-lcWGS-admixture-PAN-17vii23.csv", stringsAsFactors = F) %>% filter(K %in% 2) %>% 
#   select(1,3,4) %>% rename("ID"="Sample_ID") %>% 
#   mutate(Membership = ifelse(Q2>0.5, "BEL-WHB","HS"))  # try to put 0.6 since with 0.5 6 PAN ddRAD belugas are identified here as BEL-WHB

# Fst
wgs.fst <- read.table("./02_Results/01_ddRAD_Bringloe/03_FST/MOBELS_S_20_00703_Beluga_invariant.variant_16v23_final_average_fst.txt", header = T)


# Merge ddRAD results together
d.wgs <- wgs.pop.data %>% filter(ID %in% wgs.PCA$SampleID) %>% filter(!duplicated(ID)) %>%  # remove unused ID and duplicates
  left_join(wgs.PCA, by = c("ID"="SampleID")) %>% 
  left_join(wgs.admx, by = "ID")
colnames(d.wgs)[1] <- c("ID_GQ") 

d.wgs %>% group_by(Membership) %>% summarise(N = n())
d.wgs <- d.wgs %>% mutate(Membership = gsub("PAN","HBSC",Membership),  # WHB included with PAN
                          Membership = gsub("WHB","HBSC",Membership),  # WHB included with PAN
                          Membership = gsub("JAM-BEL","JAM",Membership),
                          Membership = gsub("CBS","CS",Membership),
                          Marker = "lcWGS")# %>% 
  # mutate(Q1.pan = NA,
  #        Q2.pan = NA,
  #        Q3.pan = NA,
  #        Membership.pan = NA,
  #        Q1.jam = NA,
  #        Q2.jam = NA,
  #        Membership.jam = NA,
  #        Membership = Membership.all)


## ddRAD + lcWGS together -------------------------------------------------

id.rad <- d.rad %>% pull(ID_GQ) %>% sort()
id.wgs <- d.wgs %>% pull(ID_GQ) %>% sort()
id.over <- id.rad[which(id.rad %in% id.wgs)]  # ID used for analyses in both datasets

str(d.rad)
str(d.wgs)
d.rad <- transform(d.rad, Lat = as.numeric(d.rad$Lat),
                          Lon = as.numeric(d.rad$Lon))
d <- bind_rows(d.rad, d.wgs) %>% select(1:33,36,34:35) %>% arrange(ID_GQ)
# over <- d %>% filter(ID_GQ %in% id.over)
# over %>% group_by(ID_GQ) %>% summarise(Cluster = Membership.all[1]==Membership.all[2]) %>% print(n = 100
# Membership of belugas in both ddRAD and lcWGS datasets is the same
d <- d %>% mutate(Marker = ifelse(ID_GQ %in% id.over, "both", Marker)) %>% 
  filter(!duplicated(ID_GQ)) %>% 
  mutate(Season = ifelse(is.na(Month), "Unknown",
                         ifelse(Month %in% c(7,8), "Summer",
                                ifelse(Month > 3 & Month < 7, "Spring",
                                       ifelse(Month > 8 & Month < 12, "Fall",
                                              "Winter")))),
         Season = factor(Season, levels = c("Spring","Summer","Fall","Winter","Unknown")))

d %>% group_by(Region1, Marker) %>% summarise(N = n()) %>% print(n = 50)  # to confirm Table 1



# PCA ---------------------------------------------------------------------

# Check all this stuff out
# # Upload lcWGS PCA results
# PCA.all.lcWGS <- read.csv("./02_Results/01_ddRAD_Bringloe/01_PCA/MOBELS-S_20_00703-Beluga-lcWGS-PCA-17vii23.csv", stringsAsFactors = F)
# PCA.all.lcWGS <- PCA.all.lcWGS %>% select(1:8,23:25)
# pop.data.wgs <- read.csv("./00_Data/03b_ddRAD_Bringloe/MOBELS_WGS_230717.csv", stringsAsFactors = F) %>% filter(!duplicated(Numero_unique_specimen))
# PCA.all.lcWGS <- PCA.all.lcWGS %>% left_join(pop.data.wgs %>% select(Numero_unique_specimen,Lieu_echantillonnage,Annee_echantillonnage,Mois_echantillonnage,
#                                                                      Jour_echantillonnage), by = c("SampleID"="Numero_unique_specimen")) %>% 
#   rename(ID = SampleID,
#          Region1 = Harvest.Region,
#          Location = Lieu_echantillonnage,
#          Year = Annee_echantillonnage,
#          Month = Mois_echantillonnage,
#          Day = Jour_echantillonnage) %>% 
#   mutate(Region1 = gsub("SAN", "BEL", Region1))
# 
# 
# # Descriptive stats
# ## SLE: proportion of belugas from St Lawrence identified as SLE
# SLE.rad <- PCA.all.MAF05NA05 %>% filter(Region1 %in% "SLE") %>% select(ID, score.PC1) %>% rename(PC1 = score.PC1) %>% 
#   mutate(marker = "RAD", ID = gsub("_rep", "", ID), cluster = ifelse(PC1>10, "SLE", "other"))
# SLE.wgs <- PCA.all.lcWGS %>% filter(Region1 %in% "SLE") %>% select(ID, PC1) %>% mutate(marker = "WGS", cluster = ifelse(PC1>0.1, "SLE", "other"))
# SLE <- bind_rows(SLE.rad,SLE.wgs)
# SLE.dup <- SLE %>% filter(duplicated(ID)) %>% pull(ID)  # ID present in both ddRAD and lcWGS datasets
# SLE <- SLE %>% mutate(marker = ifelse(ID %in% SLE.dup, "both", marker)) %>% filter(!duplicated(ID))
# SLE %>% group_by(cluster) %>% summarise(N = n())
# 
# ## JAM-BEL cluster: seasons
# JAM.rad <- PCA.all.MAF05NA05 %>% filter(score.PC2>5) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "RAD", ID = gsub("_rep", "", ID))
# JAM.wgs <- PCA.all.lcWGS %>% filter(PC2>0.1) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "WGS")
# JAM <- bind_rows(JAM.rad,JAM.wgs)
# JAM.dup <- JAM %>% filter(duplicated(ID)) %>% pull(ID)  # ID present in both ddRAD and lcWGS datasets
# JAM <- JAM %>% mutate(marker = ifelse(ID %in% JAM.dup, "both", marker)) %>% filter(!duplicated(ID))
# JAM %>% group_by(Season,Region1,marker) %>% summarise(N = n())
# JAM %>% group_by(Season,Region1,Location) %>% summarise(N = n())
# 
# ## SEH and LGR: sampling bias and cluster location
# SEH.rad <- PCA.all.MAF05NA05 %>% filter(Region1 %in% "SEH") %>% select(ID, Region1, Season, Location) %>% mutate(marker = "RAD", ID = gsub("_rep", "", ID))
# SEH.wgs <- PCA.all.lcWGS %>% filter(Region1 %in% "SEH") %>% select(ID, Region1, Season, Location) %>% mutate(marker = "WGS")
# SEH <- bind_rows(SEH.rad,SEH.wgs)
# SEH.dup <- SEH %>% filter(duplicated(ID)) %>% pull(ID)  # ID present in both ddRAD and lcWGS datasets
# SEH <- SEH %>% mutate(marker = ifelse(ID %in% SEH.dup, "both", marker)) %>% filter(!duplicated(ID))
# SEH %>% filter(Season %in% "Summer") %>% group_by(Location) %>% summarise(N = n()) %>% print(n = 50)
# 
# LGR.rad <- PCA.all.MAF05NA05 %>% filter(score.PC3 < -5) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "RAD", ID = gsub("_rep", "", ID))
# LGR.wgs <- PCA.all.lcWGS %>% filter(PC4>0.1) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "WGS")
# LGR <- bind_rows(LGR.rad,LGR.wgs)
# LGR.dup <- LGR %>% filter(duplicated(ID)) %>% pull(ID)  # ID present in both ddRAD and lcWGS datasets
# LGR <- LGR %>% mutate(marker = ifelse(ID %in% LGR.dup, "both", marker)) %>% filter(!duplicated(ID))
# LGR %>% filter(Season %in% "Summer") %>% group_by(Region1,Location,marker) %>% summarise(N = n()) %>% print(n = 50)
# 
# ## CS cluster: how many belugas from Cumberland Sound
# CS.rad <- PCA.all.MAF05NA05 %>% filter(score.PC6 < -4.729) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "RAD", ID = gsub("_rep", "", ID))
# CS.wgs <- PCA.all.lcWGS %>% filter(PC3 < -0.1) %>% select(ID, Region1, Season, Location) %>% mutate(marker = "WGS")
# CS <- bind_rows(CS.rad,CS.wgs)
# CS.dup <- CS %>% filter(duplicated(ID)) %>% pull(ID)  # ID present in both ddRAD and lcWGS datasets
# CS <- CS %>% mutate(marker = ifelse(ID %in% CS.dup, "both", marker)) %>% filter(!duplicated(ID))
# CS %>% group_by(Region1,marker,Season) %>% summarise(N = n()) %>% print(n = 50)



# Admixture ---------------------------------------------------------------

d %>% group_by(Season,Region1,Membership,Marker) %>% summarise(N = n()) %>% print(n = 150)
d %>% group_by(Marker,Membership,Region1) %>% summarise(N = n()) %>% print(n = 150)

# JAM cluster
d %>% filter(Membership %in% "JAM") %>% nrow()
d %>% filter(Membership %in% "JAM", Season %in% "Summer") %>% nrow()  # number of JAM-BEL belugas harvested/tagged during summer

## Where and when JAM belugas are harvested
d %>% filter(Membership %in% "JAM") %>% group_by(Season, Region1) %>% summarise(N = n())

## Proportion of JAM belugas harvested in Belcher Islands
d %>% filter(Region1 %in% "BEL") %>% group_by(Season, Membership) %>% summarise(N = n())

## Sex by season for JAM-BEL belugas harvested in South-East Hudson
d %>% filter(Membership %in% "JAM", Region %in% "SEH") %>% group_by(Season, Sex.qPCR) %>% summarise(N = n())
d %>% filter(Membership %in% "JAM", Season %in% "Summer") %>% group_by(Sex.qPCR) %>% summarise(N = n())  # Sex of JAM-BEL belugas harvested during summer


# LGR cluster
d %>% filter(Membership %in% "LGR") %>% nrow()

## Where and when LGR belugas are harvested
d %>% filter(Membership %in% "LGR") %>% group_by(Season, Region1) %>% summarise(N = n())


# CS cluster
d %>% filter(Membership %in% "CS") %>% nrow()
d %>% filter(Region %in% "CSB", Season %in% "Summer") %>% select(1,7:16,18,35) %>% View()

## Where and when LGR belugas are harvested
d %>% filter(Membership %in% "CS") %>% group_by(Season, Region1) %>% summarise(N = n())


# SLE cluster
d %>% filter(Membership %in% "SLE") %>% nrow()

## Location of beluga carcassess in the St Lawrence and their cluster
d %>% filter(Region1 %in% "SLE") %>% group_by(Membership,Location) %>% summarise(N = n()) %>% print(n = 50)


# # JAM-BEL cluster
# ## Which IDs are in common from both datasets?
# admx.pan <- rad.admx.pan[which(rad.admx.pan$ID_GQ %in% wgs.admx.pan$ID),]
# admx.pan <- admx.pan %>% left_join(wgs.admx.pan, by = c("ID_GQ" = "ID"))
# colnames(admx.pan)[c(2:8)] <- c("Q1.rad","Q2.rad","Q3.rad","Membership.rad","Q1.wgs","Q2.wgs","Membership.wgs")
# admx.pan <- admx.pan %>% left_join(rad.pop.data %>% filter(!duplicated(ID_GQ)))
# admx.pan %>% group_by(Membership.rad,Membership.wgs) %>% summarise(N = n())



## Sex by cluster ---------------------------------------------------------

d <- d %>% mutate(Sex.qPCR = gsub("U", NA, Sex.qPCR))
d %>% group_by(Sex.qPCR) %>% summarise(N = n())  # F = 375, M = 498, NA = 32
chisq.test(table(d$Sex.qPCR))  # X-squared = 17.33, df = 1, p-value = 3.142e-05
d %>% filter(!is.na(Sex.qPCR)) %>% group_by(Community %in% "DFO_MPO") %>% summarise(N = n())  # 32 belugas tagged by DFO, 841 harvested by Inuit communities
d %>% filter(Community %nin% "DFO_MPO") %>% group_by(Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)  # F = 362, M = 479; NA = 29
chisq.test(table(d$Sex.qPCR[d$Community %nin% "DFO_MPO"]))  # X-squared = 16.277, df = 1, p-value = 5.472e-05

d %>% filter(Community %nin% "DFO_MPO") %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)
# Membership.all Sex.qPCR     N
# CSB            F           17
# CSB            M           22
# JAM            F            9
# JAM            M           43
# JAM            NA           1
# LGR            F           33
# LGR            M           29
# LGR            NA           1
# PAN            F          280
# PAN            M          369
# PAN            NA          25
# RES            F            3
# RES            M            6
# SLE            F           20
# SLE            M           10
# SLE            NA           2

l.cluster <- split(d, d$Membership)
lapply(l.cluster, function(i){
  data <- i[i$Community %nin% "DFO_MPO",]
  t <- table(data$Sex.qPCR)
  print(chisq.test(t))
})  # RES and WHB low N (<10)






# # Admixture: FROM PAN script
# ## Identify genetic clusters ----------------------------------------------
# 
# Q.all <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_ALL_k2-10.csv", stringsAsFactors = F)
# 
# Q.res <- Q.all %>% filter(K %in% c(2:6)) %>% 
#   left_join(Q.pan, by = c("ID_GQ","K")) %>%
#   rename(Q1.all = "Q1.x", Q2.all = "Q2.x", Q3.all = "Q3.x", Q4.all = "Q4.x", Q5.all = "Q5", Q6.all = "Q6",
#          Q1.pan = "Q1.y", Q2.pan = "Q2.y", Q3.pan = "Q3.y", Q4.pan = "Q4.y") %>% 
#   select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Community, Haplo, K, Q1.all, Q2.all, Q3.all, Q4.all, Q5.all, Q6.all, Q1.pan, Q2.pan, Q3.pan, Q4.pan) %>% 
#   arrange(K, ID_GQ)
# 
# #### Q thresholds ---------------------------------------------------------
# 
# # K = 5
# hist(Q.res$Q1.all[Q.res$K %in% 5], breaks = 20)  # SLE very differentiated (Q > 0.9)
# 
# hist(Q.res$Q3.all[Q.res$K %in% 5], breaks = 20)  # CSB
# table(Q.res$Region1[Q.res$Q3.all > 0.7 & Q.res$K %in% 5])
# Q.res %>% subset(K %in% 5 & Region1 %in% "CSB") %>% View()
# Q.res %>% subset(K %in% 5) %>% View()  # 20 CSB, 1 WHB, 1 FRB (consistent with ResDoc), possible hybrids from CSB and WHB
# 
# hist(Q.res$Q4.all[Q.res$K %in% 5], breaks = 20)  # JAM very differentiated
# Q.res %>% subset(K %in% 5) %>% View()  # Q4 > 0.75 for clear JAM identification
# Q.res %>% subset(K %in% 5 & Region1 %in% "JAM") %>% View()  # All JAM Q4 > 0.9 and LON Q4 = 0.83
# Q.res %>% subset(K %in% 5 & Region1 %in% "BEL") %>% View()  # 4 SAN Q4 > 0.9 (+ 1 at 0.89). Overall gradient from 0.75 to 0.95. Difficult to distinguish from JAM
# 
# hist(Q.res$Q5.all[Q.res$K %in% 5], breaks = 20)  # LGR very differentiated (and 2 hybrids Q5 = 0.36-37)
# table(Q.res$Region1[Q.res$Q5 > 0.70 & Q.res$K %in% 5])  # 34 SEH, 11 NEH, 1 NHS, 11 SHS (0.74: 34 SEH, 10 NEH, 1 NHS, 7 SHS, no WHB; 0.6: 34 SEH, 11 NEH, 1 NHS, 12 SHS, 1 NWH)
# Q.res %>% subset(K %in% 5) %>% View()  # A few individuals from NEH, SHS, and WHB betwee 0.6 and 0.74 ( I consider those above 0.74 as 'real' LGR individuals harvested during migrations)
# Q.res %>% subset(K %in% 5 & Region1 %in% "SEH") %>% View()  #Q5 > 0.74 for clear LGR identification.
# 
# hist(Q.res$Q2.all[Q.res$K %in% 5], breaks = 20)  # PAN
# Q.res %>% subset(K %in% 3 & !is.na(Q1.pan)) %>% View()  # Q1.pan = PAN, Q2.pan = SEH, Q3.pan = WHB/NWH
# hist(Q.res$Q1.pan[Q.res$K %in% 3], breaks = 20)  # PAN.all
# hist(Q.res$Q2.pan[Q.res$K %in% 3], breaks = 20)  # SEH
# table(Q.res$Region1[Q.res$Q2.pan > 0.6 & Q.res$K %in% 3])  # 17 SEH, 1 SHS (0.5: 20 SEH, 1 SHS)
# hist(Q.res$Q3.pan[Q.res$K %in% 3], breaks = 20)  # PAN.west
# table(Q.res$Region1[Q.res$Q3.pan > 0.6 & Q.res$K %in% 3])  # 18 BEL, 3 JAM, 1 NEH, 11 NHB, 1 NHS, 27 NWH, 2 SEH, 6 SHS, 2 SWH, 1 UNG
# # (0.5: 22 BEL, 3 JAM, 1 NEH, 7 NHB, 2 NHS, 40 NWH, 9 REB, 5 SEH, 12 SHS, 3 SWH, 5 UNG)
# Q.res %>% subset(K %in% 3 & Q3.pan > 0.6) %>% group_by(Region1, Month) %>% summarise(N = n()) %>% print(n = 50)
# # Q3.pan: found at higher proportion in Western HB
# Q.res %>% subset(K %in% 3 & Q1.pan > 0.6) %>% group_by(Region1, Month) %>% summarise(N = n()) %>% print(n = 50)
# # Q1.pan: found everywhere is summertime
# Q.res %>% subset(K %in% 3 & Q3.pan < 0.5 & Q2.pan < 0.5) %>% group_by(Region1, Month) %>% summarise(N = n()) %>% print(n = 50)
# Q.res %>% subset(!is.na(Q1.pan) & K %in% 3) %>% group_by(Region1) %>% summarise(Mean_Q1 = mean(Q1.pan),
#                                                                                 Mean_Q2 = mean(Q2.pan),
#                                                                                 Mean_Q3 = mean(Q3.pan))
# # Region1 Mean_Q1 Mean_Q2 Mean_Q3
# # BEL       0.264  0.0940   0.642
# # CSB       0.726  0.0950   0.179
# # FRB       0.716  0.122    0.162
# # JAM       0.315  0.0720   0.613
# # NEH       0.574  0.176    0.250
# # NHB       0.446  0.0668   0.487
# # NHS       0.697  0.109    0.194
# # NWH       0.348  0.0845   0.568
# # SEH       0.367  0.387    0.247
# # SHS       0.602  0.128    0.270
# # SLE       0.629  0.163    0.208
# # SWH       0.598  0.108    0.294
# # UNG       0.616  0.119    0.265
# 
# Q.res.cluster.all <- Q.res %>% subset(K %in% 5) %>% 
#   mutate(Membership.all = ifelse(Q1.all > 0.7, "SLE",
#                                  ifelse(Q3.all > 0.7, "CSB",
#                                         ifelse(Q4.all > 0.7, "JAM",
#                                                ifelse(Q5.all > 0.7, "LGR",
#                                                       "PAN"))))) %>% 
#   select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Haplo, Community, K, Q1.all, Q2.all, Q3.all, Q4.all, Q5.all, Membership.all)
# 
# Q.res.cluster.pan <- Q.res %>% subset(K %in% 3 & !is.na(Q1.pan)) %>% 
#   mutate(Membership.pan = ifelse(Q2.pan > 0.5, "SEH",
#                                  ifelse(Q3.pan > 0.5, "WHB.BEL",
#                                         "PAN"))) %>% 
#   select(ID_GQ, ID_rcpt, Region1, Location, Year, Month, Day, Haplo, Community, K, Q1.pan, Q2.pan, Q3.pan, Membership.pan)
# 
# pop.admx.Q <- Q.res.cluster.all %>% 
#   left_join(Q.res.cluster.pan %>% select(ID_GQ, K, Q1.pan, Q2.pan, Q3.pan, Membership.pan), by = "ID_GQ") %>% 
#   rename(K.all = "K.x", K.pan = "K.y") %>% 
#   mutate(Membership = ifelse(is.na(Membership.pan), Membership.all, Membership.pan)) %>%  # to remove unused Pop levels
#   mutate(Season = ifelse(is.na(Month), "Unknown",
#                          ifelse(Month %in% c(7,8), "Summer",
#                                 ifelse(Month > 3 & Month < 7, "Spring",
#                                        ifelse(Month > 8 & Month < 12, "Fall",
#                                               "Winter")))))
# 
# write.csv(pop.admx.Q, file = "./02_Results/01_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture.csv", row.names = F)
# 
# 
# # Upload lcWGS admixture results + mtDNA haplo dataset
# pop.admx.Q.lcWGS <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/MOBELS-S_20_00703-Beluga-lcWGS-Metadata+K6+Het-9v23.csv", stringsAsFactors = F)
# lcWGS.meta <- read.csv("../ACCESS/MOBELS_WGS_230717.csv", stringsAsFactors = F) %>% filter(!duplicated(Numero_unique_specimen)) %>% 
#   mutate(Region_echantillonnage = gsub("SAN","BEL", Region_echantillonnage))
# 
# pop.admx.Q.lcWGS <- pop.admx.Q.lcWGS %>% select(ID,mitogenome.accession,Overlap.Luca,Membership,Q1,Q2,Q3,Q4,Q5,Q6,O.HOM.,E.HOM.,N_SITES,'F') %>% 
#   mutate(Membership = gsub("JAM-BEL", "JAM", Membership),
#          Membership = gsub("CBS", "CSB", Membership)) %>% 
#   left_join(lcWGS.meta %>% select(Numero_unique_specimen,Numero_reception_specimen,Genre,Age,Sexe_visuel,Sexe_laboratoire,Region_echantillonnage,Lieu_echantillonnage,
#                                   Latitude_echantillonnage_DD,Longitude_echantillonnage_DD,Annee_echantillonnage,Mois_echantillonnage,Jour_echantillonnage,
#                                   Communaute,Haplotype), by = c("ID"="Numero_unique_specimen")) %>% 
#   # select(ID,Region,Location,Community,Year,Month,Day,Season,Haplo,Membership) %>% 
#   # rename("ID_GQ"="ID","Region1" = "Region","Membership.all"="Membership") %>%
#   mutate(Marker = "lcWGS",
#          Membership.all = Membership) %>%  # since ddRAD dataset has both columns
#   mutate(Season = ifelse(is.na(Mois_echantillonnage), "Unknown",
#                          ifelse(Mois_echantillonnage %in% c(7,8), "Summer",
#                                 ifelse(Mois_echantillonnage > 3 & Mois_echantillonnage < 7, "Spring",
#                                        ifelse(Mois_echantillonnage > 8 & Mois_echantillonnage < 12, "Fall",
#                                               "Winter")))))
# colnames(pop.admx.Q.lcWGS)[c(1,15:28)] <- c("ID_GQ","ID_rcpt","Species","Age","Sex.visual","Sex.qPCR","Region1","Location","Lat","Lon","Year","Month","Day","Community","Haplo") 
# 
# # Merge lcWGS and ddRAD admixture results
# pop.admx.all <- pop.admx.Q %>% select(ID_GQ,Region1,Location,Community,Year,Month,Day,Season,Haplo,Membership.all,Membership) %>% 
#   mutate(ID_GQ = gsub("_rep", "", ID_GQ),
#          Marker = "ddRAD") %>% 
#   bind_rows(pop.admx.Q.lcWGS %>% select(ID_GQ,Region1,Location,Community,Year,Month,Day,Season,Haplo,Membership.all,Membership,Marker))
# dup.id <- pop.admx.all %>% filter(duplicated(ID_GQ)) %>% pull(ID_GQ)  # duplicated ID (73 IDs)
# # pop.admx.all %>% filter(ID_GQ %in% dup.id) %>% group_by(ID_GQ,Membership.all) %>% summarise(N = n()) %>% print(n = 150)  # all duplicates have same admixture membership
# pop.admx.all <- pop.admx.all %>% arrange(ID_GQ) %>%  mutate(Marker = ifelse(ID_GQ %in% dup.id, "both", Marker)) %>% filter(!duplicated(ID_GQ))  # specify marker for duplicated IDs and remove duplicates
# pop.admx.all %>% group_by(Marker) %>% summarise(N = n())
# pop.admx.all %>% group_by(Region1, Marker) %>% summarise(N = n()) %>% print(n = 50)
# pop.admx.all %>% filter(Region1 %nin% c("RES","SLE")) %>% group_by(Membership.all) %>% summarise(N = n())
# pop.admx.all %>% group_by(Membership.all) %>% summarise(N = n())  # 677 out of 862 belugas are PAN
# pop.admx.all %>% group_by(Membership) %>% summarise(N = n())
# 
# 

# 
# 
# 
# 
# # Misc descriptive stats: LGR cluster temporal harvest
# pop.admx.all %>% filter(Membership %in% "LGR") %>% group_by(Region1, Location, Year, Month) %>% summarise(N = n()) %>% print(n = 50)
# 
# 
# # Misc descriptive stats for association admixture membership & mtDNA haplo for LGR, SEH, and PAN harvested in SEH during summer
# h.ehb <- c("HL041","HL049","HL112","HL002","HL097","HL071","HL059","HL043","HL005","HL031","HL007","HL069","HL017","HL020","HL083","HL034","HL038","HL018",
#            "HL063","HL019","HL016","HL046","HL122","HL033","HL099","HL102","HL039","HL129","HL077","HL103","HL092","HL10","HL009","HL123")  # HL099 + HL102 SLE
# pop.admx.all %>% filter(Membership %in% "LGR") %>% group_by(Region1, Season, Haplo %in% h.ehb) %>% summarise(N = n()) %>% print(n = 100)  # 11 haplo from western and 19 from eastern haplogroup
# pop.admx.all %>% filter(Membership %in% "LGR", Season %in% "Summer", Region1 %in% "SEH") %>% group_by(Haplo) %>% summarise(N = n())
# pop.admx.all %>% filter(Membership %in% "SEH") %>% group_by(Region1, Season, Haplo %in% h.ehb) %>% summarise(N = n()) %>% print(n = 100)  # 6 haplo from western and 11 from eastern haplogroup
# pop.admx.all %>% filter(Membership %in% "SEH", Season %in% "Summer", Region1 %in% "SEH") %>% group_by(Haplo) %>% summarise(N = n())
# # Reminder for PAN: this includes possible BEL.WHB from lcWGS dataset
# pop.admx.all %>% filter(Membership %in% "PAN") %>% group_by(Region1, Season, Haplo %in% h.ehb) %>% summarise(N = n()) %>% print(n = 100)  # 18 haplo from western and 40 from eastern haplogroup
# pop.admx.all %>% filter(Membership %in% "PAN", Season %in% "Summer", Region1 %in% "SEH") %>% group_by(Haplo) %>% summarise(N = n())
# 
# # Tagged belugas (LGR)
# pop.admx.all %>% filter(Community %in% "DFO_MPO") %>% group_by(Season, Region1, Location, Membership) %>% summarise(N = n())
# pop.admx.all %>% filter(Community %in% "DFO_MPO") %>% pull(ID_GQ)  # tagged belugas
# # "S_20_00600" "S_20_00601" "S_20_00602" "S_20_00620" "S_20_00621" "S_20_00622" "S_20_00652" "S_20_00653" "S_20_00699" "S_20_00710" "S_20_00711" "S_20_00712"
# # "S_20_00901" "S_20_00903" "S_20_00904" "S_20_01021" "S_20_02173" "S_20_02174" "S_20_02175" "S_20_02177" "S_20_02178" "S_20_02179" "S_20_02180" "S_20_03444"
# # "S_20_03446" "S_20_03447" "S_20_03448" "S_20_03449" "S_20_03450" "S_20_03451" "S_20_03452" "S_20_03453" "S_20_03454" "S_20_03455" "S_20_03507"
# pop.admx.all %>% filter(Community %in% "DFO_MPO", Membership %in% "LGR") %>% pull(ID_GQ)  # tagged LGR belugas
# # "S_20_00600" "S_20_00601" "S_20_00602" "S_20_00620" "S_20_00621" "S_20_00901" "S_20_01021" "S_20_02174" "S_20_02175" "S_20_02177" "S_20_02178" "S_20_02179"
# # "S_20_02180"
# # 9 of the 17 EHB in Bailleul et al. 2012 followed till wintertime
# 
# # Misc descriptive stats for association cluster - season - harvest region/location
# pop.admx.all %>% filter(Membership %in% "LGR") %>% group_by(Season, Month, Year, Day, Region1, Location) %>% summarise(N = n()) %>% print(n = 100)
# lgr.mig <- pop.admx.all %>% filter(Membership %in% "LGR", Season %nin% c("Summer","Unknown")) %>% select(ID_GQ, Region1, Location, Season, Marker) %>% arrange(desc(Season), Region1)
# lgr.mig %>% group_by(Season, Region1, Location) %>% summarise(N = n())
# 
# pop.admx.all %>% filter(Membership %in% "JAM") %>% group_by(Season, Month, Year, Day, Region1, Location) %>% summarise(N = n()) %>% print(n = 100)
# 
# pop.admx.all %>% filter(Membership %in% "CSB") %>% group_by(Season, Month, Year, Day, Region1, Location) %>% summarise(N = n()) %>% print(n = 100)
# 
# pop.admx.all %>% filter(Membership.all %in% "PAN") %>% group_by(Season, Region1, Location) %>% summarise(N = n()) %>% print(n = 200)  # includes BEL.WHB and SEH
# 
# pop.admx.all %>% filter(Membership %in% "WHB.BEL") %>% group_by(Season,Month, Region1, Location) %>% summarise(N = n()) %>% print(n = 100)
# 
# pop.admx.all %>% filter(Membership %in% "SEH") %>% group_by(Season,Month, Region1, Location) %>% summarise(N = n()) %>% print(n = 100)
# 
# pop.admx.Q %>% filter(Membership %in% "PAN") %>% group_by(Season, Month, Region1, Location) %>% summarise(N = n()) %>% print(n = 200)  # includes BEL.WHB from lcWGS dataset
# 





# FST ---------------------------------------------------------------------

head(rad.fst)
rad.fst <- rad.fst %>% mutate(Population1 = gsub("CSB","CS",Population1),
                              Population2 = gsub("CSB","CS",Population2))
head(wgs.fst)
names(wgs.fst) <- c("Population1", "Population2", "Fst")
unique(wgs.fst$Population1)  # remove PAN.west (IDs chosen according to ddRAD results) + BEL (described according lcWGS Admixture results)
wgs.fst <- wgs.fst %>% mutate(Population1 = gsub("CBS","CS",Population1),
                              Population2 = gsub("CBS","CS",Population2),
                              Population1 = gsub("PAN.all","HBSC",Population1),
                              Population2 = gsub("PAN.all","HBSC",Population2)) %>% 
  arrange(Population1, Population2) %>% 
  filter(!Population1 %in% c("BEL","PAN.west")) %>% 
  filter(!Population2 %in% c("BEL","PAN.west"))

gFST.lcWGS.6Q <- wgs.fst %>%
  mutate(Population1 = factor(Population1, levels = c("RES","HBSC","JAM","LGR","CS","SLE")),
         Population2 = factor(Population2, levels = c("RES","HBSC","JAM","LGR","CS","SLE"))) %>%
  ggplot(aes(x=Population1, y=Population2, fill=Fst)) +
  geom_tile(colour=  "white") +
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
gFST.lcWGS.6Q

## Mean + ranges lcWGS -----------------------------------------------------------

pop <- sort(unique(wgs.fst$Population1))

for(i in pop){
  print(wgs.fst %>% filter(Population1 %in% i) %>% summarise(Mean = mean(Fst),
                                                           Range = paste(min(Fst), max(Fst), sep = "-"),
                                                           Pop = i))
}
#       Mean                   Range Pop
# 0.03352584 0.009911285-0.080455231  CS
#       Mean                 Range  Pop
# 0.02811875 0.009690069-0.0731839 HBSC
#       Mean                 Range Pop
# 0.05328032 0.032936107-0.1030176 JAM
#       Mean                  Range Pop
# 0.03388291 0.009690069-0.08029682 LGR
#       Mean                   Range Pop
# 0.03819738 0.014872376-0.088529121 RES
#       Mean               Range Pop
# 0.08509653 0.0731839-0.1030176 SLE




## Figures ----------------------------------------------------------------

FST.6Q <- wgs.fst %>% 
  left_join(rad.fst, by = c("Population1", "Population2")) %>%
  rename(c(Fst.lcWGS.6Q = "Fst.x", Fst.ddRAD.5Q = "Fst.y")) %>%  # , Fst.ddRAD.red.5Q = "Fst"
  mutate(Pairwise = paste0(Population1, " vs ", Population2)) %>%
  mutate(Fst.6Q = ifelse(Pairwise %in% "RES vs HBSC", Fst.lcWGS.6Q[Pairwise %in% "RES vs HBSC"],
                  ifelse(Pairwise %in% "HBSC vs RES", Fst.ddRAD.5Q[Pairwise %in% "HBSC vs RES"],
                  ifelse(Pairwise %in% "RES vs JAM", Fst.lcWGS.6Q[Pairwise %in% "RES vs JAM"],
                  ifelse(Pairwise %in% "JAM vs RES", Fst.ddRAD.5Q[Pairwise %in% "JAM vs RES"],
                  ifelse(Pairwise %in% "RES vs LGR", Fst.lcWGS.6Q[Pairwise %in% "RES vs LGR"],
                  ifelse(Pairwise %in% "LGR vs RES", Fst.ddRAD.5Q[Pairwise %in% "LGR vs RES"],
                  ifelse(Pairwise %in% "RES vs CS", Fst.lcWGS.6Q[Pairwise %in% "RES vs CS"],
                  ifelse(Pairwise %in% "CS vs RES", Fst.ddRAD.5Q[Pairwise %in% "CS vs RES"],
                  ifelse(Pairwise %in% "RES vs SLE", Fst.lcWGS.6Q[Pairwise %in% "RES vs SLE"],
                  ifelse(Pairwise %in% "SLE vs RES", Fst.ddRAD.5Q[Pairwise %in% "SLE vs RES"],
                         ifelse(Pairwise %in% "HBSC vs JAM", Fst.lcWGS.6Q[Pairwise %in% "HBSC vs JAM"],
                         ifelse(Pairwise %in% "JAM vs HBSC", Fst.ddRAD.5Q[Pairwise %in% "JAM vs HBSC"],
                         ifelse(Pairwise %in% "HBSC vs LGR", Fst.lcWGS.6Q[Pairwise %in% "HBSC vs LGR"],
                         ifelse(Pairwise %in% "LGR vs HBSC", Fst.ddRAD.5Q[Pairwise %in% "LGR vs HBSC"],
                         ifelse(Pairwise %in% "HBSC vs CS", Fst.lcWGS.6Q[Pairwise %in% "HBSC vs CS"],
                         ifelse(Pairwise %in% "CS vs HBSC", Fst.ddRAD.5Q[Pairwise %in% "CS vs HBSC"],
                         ifelse(Pairwise %in% "HBSC vs SLE", Fst.lcWGS.6Q[Pairwise %in% "HBSC vs SLE"],
                         ifelse(Pairwise %in% "SLE vs HBSC", Fst.ddRAD.5Q[Pairwise %in% "SLE vs HBSC"],
                                ifelse(Pairwise %in% "JAM vs LGR", Fst.lcWGS.6Q[Pairwise %in% "JAM vs LGR"],
                                ifelse(Pairwise %in% "LGR vs JAM", Fst.ddRAD.5Q[Pairwise %in% "LGR vs JAM"],
                                ifelse(Pairwise %in% "JAM vs CS", Fst.lcWGS.6Q[Pairwise %in% "JAM vs CS"],
                                ifelse(Pairwise %in% "CS vs JAM", Fst.ddRAD.5Q[Pairwise %in% "CS vs JAM"],
                                ifelse(Pairwise %in% "JAM vs SLE", Fst.lcWGS.6Q[Pairwise %in% "JAM vs SLE"],
                                ifelse(Pairwise %in% "SLE vs JAM", Fst.ddRAD.5Q[Pairwise %in% "SLE vs JAM"],
                                       ifelse(Pairwise %in% "LGR vs CS", Fst.lcWGS.6Q[Pairwise %in% "LGR vs CS"],
                                       ifelse(Pairwise %in% "CS vs LGR", Fst.ddRAD.5Q[Pairwise %in% "CS vs LGR"],
                                       ifelse(Pairwise %in% "LGR vs SLE", Fst.lcWGS.6Q[Pairwise %in% "LGR vs SLE"],
                                       ifelse(Pairwise %in% "SLE vs LGR", Fst.ddRAD.5Q[Pairwise %in% "SLE vs LGR"],
                                              ifelse(Pairwise %in% "CS vs SLE", Fst.lcWGS.6Q[Pairwise %in% "CS vs SLE"],
                                              ifelse(Pairwise %in% "SLE vs CS", Fst.ddRAD.5Q[Pairwise %in% "SLE vs CS"],
                                                     NA)))))))))))))))))))))))))))))))

### Matrix ----------------------------------------------------------------

gFST.6Q <-  FST.6Q %>% filter(!is.na(Fst.6Q)) %>%
  mutate(Population1 = factor(Population1, levels = c("RES","HBSC","JAM","LGR","CS","SLE")),
         Population2 = factor(Population2, levels = c("RES","HBSC","JAM","LGR","CS","SLE"))) %>%
  ggplot(aes(x=Population1, y=Population2, fill=Fst.6Q)) +
  geom_tile(colour=  "white") +
  geom_text(aes(label = round(Fst.6Q, digits = 3)), color = "white", size = 9) +  # fontface = "bold"
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
        axis.text.x = element_text(angle =  0, vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 27.5, color = "black"))
pdf(file = "02_Results/01_ddRAD_Bringloe/03_FST/FST.6Q.lcWGS.ddRAD.240212.pdf", width = 11.5, height = 9.25)
gFST.6Q
dev.off()


### Correlation -----------------------------------------------------------

FST.5Q.cor <- FST.6Q %>%  # 5Q because the correlation is between 5 clusters (6th is RES not present in ddRAD dataset)
  filter(!duplicated(str_c(pmin(Population1, Population2), pmax(Population1,Population2))))  # from https://stackoverflow.com/questions/32035865/r-remove-rows-from-a-data-frame-that-contain-a-duplicate-of-either-combination-o

cor.test(FST.5Q.cor$Fst.lcWGS.6Q, FST.5Q.cor$Fst.ddRAD.5Q, method = "s")  # rho = 0.988, p < 0.001

gFST.cor.5Q <- FST.5Q.cor %>% ggplot(aes(x = Fst.ddRAD.5Q, y = Fst.lcWGS.6Q, colour = Pairwise)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, col = "black", se = F) +
  geom_point(size = 14, stroke = 1, col = alpha("black", 0.6), fill = alpha("grey", 0.6), shape = 21) +
  scale_x_continuous(limits = c(0,0.135), breaks = c(0,0.025,0.05,0.075,0.1,0.125,0.15)) +
  scale_y_continuous(limits = c(0,0.135), breaks = c(0,0.025,0.05,0.075,0.1,0.125)) +
  labs(x = "Fst ddRAD", y = "Fst lcWGS") +
  theme_bw() +
  theme(axis.text = element_text(size = 25, colour = "black"),
        axis.title.y.right = element_text(angle = 90),
        axis.title = element_text(size = 30, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
gFST.cor.5Q <- gFST.cor.5Q + annotate("text", label = substitute(paste(italic("r"), " = 0.988")), x = 0.056, y = 0.061, size = 10, angle = 48)
# pdf(file = "02_Results/01_ddRAD_Bringloe/03_FST/FST.5Q.cor.lcWGS.ddRAD.240212.pdf", width = 10, height = 9)
gFST.cor.5Q
dev.off()


