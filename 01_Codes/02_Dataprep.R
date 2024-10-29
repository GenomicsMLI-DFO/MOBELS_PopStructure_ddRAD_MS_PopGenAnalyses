# Info --------------------------------------------------------------------
#
# Overview: Prepare dataset for eastern Arctic beluga population structure
# 
# Author: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2022-10-19
#


# Housekeeping ------------------------------------------------------------

# Verify if you're in the right directory
getwd()

# Clear workspace
rm(list = ls())

# Libraries
# if(!require(tidyverse)){install.packages("tidyverse")}
library(readxl)
library(dplyr)
# library(ggplot2)
# library(tidyr)

# Functions
"%nin%" <- Negate("%in%")



# Data --------------------------------------------------------------------

## IBIS
IBIS.barcode <- read_excel("../Beluga_ddRAD_Novaseq_2019/00_Data/00_FileInfos/BarcodesIBIS.xlsx")
IBIS.barcode

## GQ
d <- read.csv("../ACCESS/MOBELS_ddRAD811_230717.csv", stringsAsFactors = F)
d <- d %>% mutate(Cat_sample = ifelse(duplicated(Numero_unique_specimen ), "Duplicate", "Sample"),
                  ID_GQ = ifelse(Cat_sample == "Duplicate", paste0(Numero_unique_specimen, "_rep"), Numero_unique_specimen)) %>% 
  mutate(ID_GQ = ifelse(duplicated(ID_GQ), paste0(ID_GQ, "_1"), ID_GQ)) %>% 
  arrange(ID_GQ)


## Merge datasets ---------------------------------------------------------

# Merge d and GQ.data: dataset for ddRAD analyses
rad <- d %>% 
  left_join(IBIS.barcode %>% select(No_puits_envoi = Cell2, Barcode))




# Data inspection ---------------------------------------------------------

## ddRAD ------------------------------------------------------------------

rad %>% View()
rad <- rad %>% select(ID_GQ,Numero_unique_specimen,Numero_unique_extrait,Numero_reception_specimen,Cat_sample,No_plaque_envoi,No_puits_envoi,Barcode,
                      Nom_commun,Genre,Espece,Age,Sexe_visuel,Sexe_laboratoire,Region_echantillonnage,Lieu_echantillonnage,Latitude_echantillonnage_DD,
                      Longitude_echantillonnage_DD,Annee_echantillonnage,Mois_echantillonnage,Jour_echantillonnage,Communaute,Haplotype)
colnames(rad) <- c("ID_GQ","ID","ID_Ext","ID_rcpt","Cat_sample","No_plate","No_well","Barcode","Common_name","Genus","Species","Age","Sex_visual","Sex_qPCR",
                   "Region","Location","Lat","Lon","Year","Month","Day","Community","Haplo")

rad %>% filter(Cat_sample %in% "Sample") %>%
  group_by(Region) %>%
  summarise(N = n()) %>% View()

table(rad$Region, useNA = "ifany")  # previously NHB and REB not otgether + 6 IDs switched from UNG to SHS
# CSB FRB JAM LON NEH NHB NHS NWH SAN SEH SHS SLE SWH UNG 
#  42  20  25   1  40  58  30  80  54 144 144  33  15 125

rad <- rad %>%
  mutate(Region1 = ifelse(Region %in% "CSB", "CSB",
                          ifelse(Region %in% "FRB", "FRB",
                                 ifelse(Region %in% c("JAM","LON"), "JAM",
                                        ifelse(Region %in% "NEH", "NEH",
                                               ifelse(Region %in% "NHB", "NHB",
                                                      ifelse(Region %in% "NHS", "NHS",
                                                             ifelse(Region %in% "SAN", "BEL",
                                                                    ifelse(Region %in% "SEH", "SEH",
                                                                           ifelse(Region %in% "SHS", "SHS",
                                                                                  ifelse(Region %in% "SLE", "SLE",
                                                                                         ifelse(Region %in% "UNG", "UNG",
                                                                                                ifelse(Region %in% "NWH", "NWH",
                                                                                                       ifelse(Region %in% "SWH", "SWH",
                                                                                                                     NA)))))))))))))) %>% 
  mutate(Region1 = factor(Region1, levels = c("SLE","CSB","FRB","UNG","NHS","SHS","NEH","SEH","BEL","JAM","SWH","NWH","NHB"))) %>% 
  mutate(Region2 = ifelse(Region %in% "CSB", "CSB",
                          ifelse(Region %in% "FRB", "FRB",
                                 ifelse(Region %in% c("JAM","LON"), "JAM",
                                        ifelse(Region %in% "NEH", "NEH",
                                               ifelse(Region %in% c("NHB","NWH","SWH"), "WHB",
                                                      ifelse(Region %in% "NHS", "NHS",
                                                             ifelse(Region %in% "SAN", "BEL",
                                                                    ifelse(Region %in% "SEH", "SEH",
                                                                           ifelse(Region %in% "SHS", "SHS",
                                                                                  ifelse(Region %in% "SLE", "SLE",
                                                                                         ifelse(Region %in% "UNG", "UNG",
                                                                                                NA)))))))))))) %>% 
  mutate(Region2 = factor(Region2, levels = c("SLE","CSB","FRB","UNG","NHS","SHS","NEH","SEH","BEL","JAM","WHB")))

table(rad$Region1, useNA = "ifany")
# SLE CSB FRB UNG NHS SHS NEH SEH BEL JAM SWH NWH NHB 
#  33  42  20 125  30 144  40 144  54  26  15  80  58
table(rad$Region2, useNA = "ifany")
# SLE CSB FRB UNG NHS SHS NEH SEH BEL JAM WHB 
#  33  42  20 125  30 144  40 144  54  26 153

table(rad$Year, useNA = "ifany")  # all good
table(rad$Month, useNA = "ifany")  # 0 should be NA
#   5    6    7    8    9   10   11   12 <NA> 
#  36  125  267  185   50  104   12    4   28
# rad$Month[rad$Month %in% 0] <- NA




# Save datasets -----------------------------------------------------------

if(!file.exists(file.path("./00_Data", "02_Dataset"))){
  dir.create(file.path("./00_Data", "02_Dataset"))
  print(file.path("./00_Data", "02_Dataset"))
}

write.csv(rad, file = file.path("./00_Data/02_Dataset/", "Beluga_ddRAD.csv"), row.names = F)
