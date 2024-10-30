library(tidyverse)
options(scipen = 999)

`%nin%` <- Negate(`%in%`)



# Data --------------------------------------------------------------------

## Upload -----------------------------------------------------------------
data.rad <- read.csv("./00_Data/02_Dataset/MOBELS-S_20_00703-Beluga-ddRADseq-Metadata+K5_19ix23.csv")
data.wgs <- read.csv("./00_Data/02_Dataset/lcWGS/MOBELS_WGS_230717.csv")
data.wgs.mt <- read.csv("./00_Data/02_Dataset/lcWGS/MOBELS-mitogenomes-lcWGS-PCA-20ix23.csv", na.string = c("","NA"))

# meta.rad.red <- read.csv("./00_Data/03b_ddRAD_Bringloe/Beluga_ddRAD_analyses.csv", stringsAsFactors = F)
res.wgs.admx.het <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/lcWGS/MOBELS-S_20_00703-Beluga-lcWGS-Metadata+K6+Het-19ix23.csv")
# res.wgs.admx.pan <- read.csv("./02_Results/01b_ddRAD_Bringloe/02_Admixture/MOBELS-S_20_00703-Beluga-lcWGS-admixture-PAN-17vii23.csv")
# res.wgs.admx.jam <- read.csv("./02_Results/01b_ddRAD_Bringloe/02_Admixture/MOBELS-S_20_00703-Beluga-lcWGS-admixture-JAM-BEL-17vii23.csv")

res.rad.admx <- read.csv("./02_Results/01_ddRAD_Bringloe/02_Admixture/02b_neutral_SNPs/Beluga_ddRAD_meta_admixture_neutral_K5.csv")
# res.rad.admx.pan <- read.csv("./02_Results/01b_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_PAN_K3.csv")
# res.rad.admx.jam <- read.csv("./02_Results/01b_ddRAD_Bringloe/02_Admixture/Beluga_ddRAD_meta_admixture_JAM_K2.csv")
res.rad.het <- read.delim("./02_Results/01_ddRAD_Bringloe/04_Heterozygosity/04b_neutral_SNPs/out.het", skip=0, sep = "\t", header = T )


## Keep popgen samples only -----------------------------------------------

d.rad <- data.rad %>% mutate(ID_GQ = ifelse(duplicated(Numero_unique_specimen), paste0(Numero_unique_specimen, "_rep"),Numero_unique_specimen),
                             ID_GQ = ifelse(duplicated(ID_GQ), paste0(ID_GQ, "_1"), ID_GQ)) %>%  # update names of duplicated (triplicated) specimens
  filter(ID_GQ %in% res.rad.admx$ID_GQ) %>%  # keep only specimens used for popgen analyses
  select(ID_GQ, BioSample.Acession, Type_analyse, Sexe_laboratoire, Age, Region_echantillonnage, Lieu_echantillonnage, Annee_echantillonnage, Mois_echantillonnage, Jour_echantillonnage) %>% 
  mutate(Type_analyse = ifelse(gsub("_rep.*", "", ID_GQ) %in% res.wgs.admx.het$ID, "Both", Type_analyse),
         Sexe_laboratoire = ifelse(Sexe_laboratoire %in% "U", NA, Sexe_laboratoire))
colnames(d.rad) <- c("ID","BioSample.Accession","Marker","Sex.qPCR","Age","Region","Location","Year","Month","Day")
View(d.rad)
str(d.rad)

d.wgs <- res.wgs.admx.het %>% select(ID, BioSample.Accession, mitogenome.accession, Sex.qPCR, Age, Region, Location, Year, Month, Day, 
                                     Membership, Q1, Q2, Q3, Q4, Q5, Q6, O.HOM., E.HOM., 'F') %>% 
  mutate(Marker = ifelse(ID %in% gsub("_rep.*", "", d.rad$ID), "Both", "lcWGS")) %>% 
  left_join(data.wgs.mt %>% filter(!is.na(SampleID)) %>% select(SampleID,Mitochondrial.clade), by = c("ID"="SampleID"))
colnames(d.wgs)[c(2:3,11:20)] <- c("BioSample.Accession","Mitogenome.accession","Membership.wgs","Q1.wgs","Q2.wgs","Q3.wgs",
                                   "Q4.wgs","Q5.wgs","Q6.wgs","O.Hom.wgs","E.Hom.wgs","F.wgs")
View(d.wgs)
str(d.wgs)


## Include results --------------------------------------------------------

d.rad <- d.rad %>% left_join(res.rad.admx %>% select(ID_GQ,Membership.neutral,Q1,Q2,Q3,Q4,Q5), by = c("ID"="ID_GQ")) %>% 
  left_join(res.rad.het %>% select(INDV,O.HOM.,E.HOM.,'F'), by = c("ID"="INDV"))
colnames(d.rad)[c(11:19)] <- c("Membership.rad","Q1.rad","Q2.rad","Q3.rad","Q4.rad","Q5.rad",
                               "O.Hom.rad","E.Hom.rad","F.rad")
str(d.rad)


## Combine tables ---------------------------------------------------------

d <- d.rad %>% full_join(d.wgs, by = c("ID","Marker","BioSample.Accession","Sex.qPCR","Age","Region","Location","Year","Month","Day")) %>% 
  select(1,3,2,20,31,4:11,21,12:16,22:27,17:19,28:30) %>% 
  mutate(Region = gsub("SAN", "BEL", Region),
         Region = gsub("LON", "JAM", Region),
         Region = gsub("CSB", "CS", Region)) %>% 
  mutate(Membership.rad = ifelse(is.na(Membership.rad), NA,
                                 ifelse(Membership.rad %in% "PAN", "HBSC",
                                        ifelse(Membership.rad %in% "CSB", "CS", 
                                               ifelse(Membership.rad %in% "JAM", "JB",
                                                      Membership.rad)))), 
         Membership.wgs = ifelse(is.na(Membership.wgs), NA,
                                 ifelse(Membership.wgs %in% "CBS", "CS",
                                        ifelse(Membership.wgs %in% c("PAN","WHB"), "HBSC",
                                               ifelse(Membership.wgs %in% "JAM-BEL", "JB", Membership.wgs))))) %>%
  transform(Region = factor(Region, levels = c("SLE","CS","FRB","UNG","NHS","SHS","NEH","SEH","BEL","LON","JAM","SWH","NWH","NHB","RES"))) %>% 
  arrange(Region, ID)
str(d)

d <- d %>% 
  mutate(Region = ifelse(Region %in% "SLE", "Saint Lawrence",
                         ifelse(Region %in% "CS", "Cumberland Sound",
                                ifelse(Region %in% "FRB", "Frobisher Bay",
                                       ifelse(Region %in% "UNG", "Ungava Bay",
                                              ifelse(Region %in% "NHS", "North Hudson Strait",
                                                     ifelse(Region %in% "SHS", "South Hudson Strait",
                                                            ifelse(Region %in% "NEH", "North East Hudson Bay",
                                                                   ifelse(Region %in% "SEH", "South East Hudson Bay",
                                                                          ifelse(Region %in% "BEL", "Belcher Islands",
                                                                                 ifelse(Region %in% "JAM", "James Bay",
                                                                                        ifelse(Region %in% "SWH", "South West Hudson Bay",
                                                                                               ifelse(Region %in% "NWH", "North Western Hudson Bay",
                                                                                                      ifelse(Region %in% "NHB", "North Hudson Bay",
                                                                                                             ifelse(Region %in% "RES", "Resolute Bay",
                                                                                                                    NA)))))))))))))))


# write.csv(d, file = "./00_Data/02_Dataset/TableS1_241030.csv", row.names = F)

table(d$Marker)

# over <- d.rad[d.rad$Numero_unique_specimen %in% id.wgs,]
# length(which(id.rad %in% id.wgs))
# sort(id.rad[id.rad %in% id.wgs])

sex <- d %>% filter(Region %nin% "SLE", !is.na(Sex.qPCR)) %>% group_by(Sex.qPCR) %>% summarise(N = n())
chisq.test(sex$N)
# X-squared = 21.671, df = 1, p-value = 0.000003237

d %>% group_by(Month, Region) %>% summarize(N = n()) %>% print(n = 100)
# Month Region     N
# 1     BEL        1
# 2     BEL        4
# 4     SLE        1
# 5     SLE        5
# 5     NHS        2
# 5     NEH        2
# 5     BEL       42
# 6     SLE        3
# 6     CS         1
# 6     FRB        8
# 6     UNG       21
# 6     NHS        6
# 6     SHS       57
# 6     SEH       12
# 6     BEL       43
# 7     SLE       13
# 7     CS        40
# 7     FRB       11
# 7     UNG       63
# 7     NHS        7
# 7     SHS       35
# 7     SEH       51
# 7     BEL       14
# 7     JAM        4
# 7     SWH       11
# 7     NWH       17
# 7     NHB       12
# 8     SLE        4
# 8     CS         6
# 8     FRB        5
# 8     UNG        7
# 8     SEH       59
# 8     JAM       15
# 8     SWH        3
# 8     NWH       63
# 8     NHB       34
# 8     RES        6
# 9     SLE        2
# 9     UNG        1
# 9     NEH        1
# 9     SEH        4
# 9     BEL        7
# 9     JAM        5
# 9     NWH       13
# 9     NHB       15
# 9     RES        3
# 10    SLE        4
# 10    NHS       19
# 10    SHS       42
# 10    NEH       28
# 10    SEH        8
# 10    BEL        3
# 10    JAM        1
# 11    SLE        2
# 11    NHS       11
# 11    SHS       13
# 12    BEL        3
# NA    CS         1
# NA    SHS        3
# NA    NEH        1
# NA    SEH       14
# NA    BEL        8
# NA    NWH        3
# NA    NHB        2
d %>% subset(Month %in% c(7,8)) %>% group_by(Region) %>% summarise(N = n()) # summer
d %>% subset(Month %in% c(7,8)) %>% group_by(Sex.qPCR) %>% summarise(N = n()) # summer
d %>% subset(Month %in% c(3:6,9:11)) %>% group_by(Region) %>% summarise(N = n())
d %>% subset(Month %in% c(3:6,9:11)) %>% nrow()
d %>% filter(Month %in% c(7,8), Region %in% "CS") %>% group_by(Membership.rad, Membership.wgs, Sex.qPCR) %>% summarise(N = n()) # summer
d %>% filter(Region %nin% "SLE") %>% group_by(Sex.qPCR) %>% summarise(N = n())  # harvest sex ratio (but here there are tags...)

# Dispersal
d_dispersal <- d %>% filter(ifelse(Region %in% c("SLE","RES"), Month %in% c(1:12),
                                   Month %in% 7:8))  # 500 rows - Anse au Loup carcass = 499

hbsc_disp <- d_dispersal[which(d_dispersal$Membership.rad %in% "HBSC" | d_dispersal$Membership.wgs %in% "HBSC"),]
hbsc_disp %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)

cs_disp <- d_dispersal[which(d_dispersal$Membership.rad %in% "CS" | d_dispersal$Membership.wgs %in% "CS"),]
cs_disp %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)

lgr_disp <- d_dispersal[which(d_dispersal$Membership.rad %in% "LGR" | d_dispersal$Membership.wgs %in% "LGR"),]
lgr_disp %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)

jb_disp <- d_dispersal[which(d_dispersal$Membership.rad %in% "JAM" | d_dispersal$Membership.wgs %in% "JAM"),]
jb_disp %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)


# By cluster
lgr <- d[which(d$Membership.rad %in% "LGR" | d$Membership.wgs %in% "LGR"),]
lgr %>% subset(Month %in% c(7,8)) %>% nrow()
lgr %>% subset(Month %in% c(7,8)) %>% group_by(Region, Sex.qPCR) %>% summarise(N = n())
lgr %>% subset(Month %in% c(7,8)) %>% group_by(Region, Location, Sex.qPCR) %>% summarise(N = n())
lgr %>% subset(Month %in% c(4:6)) %>% group_by(Region, Month, Sex.qPCR) %>% summarise(N = n())
lgr %>% subset(Month > 8) %>% group_by(Region, Month, Sex.qPCR) %>% summarise(N = n())

cs <- d[which(d$Membership.rad %in% "CS" | d$Membership.wgs %in% "CS"),]
cs %>% subset(Month %in% c(7,8)) %>% nrow()
cs %>% subset(Month %in% c(7,8)) %>% group_by(Region, Sex.qPCR) %>% summarise(N = n())
cs %>% subset(Month %nin% c(7,8)) %>% group_by(Month, Region, Sex.qPCR) %>% summarise(N = n())

jam <- d[which(d$Membership.rad %in% "JAM" | d$Membership.wgs %in% "JAM"),]
jam %>% group_by(Month) %>% summarise(N = n())
jam %>% subset(Month %in% c(7,8)) %>% group_by(Region, Sex.qPCR) %>% summarise(N = n())
jam %>% subset(Month %in% c(9:11)) %>% group_by(Region, Month, Day, Location, Sex.qPCR) %>% summarise(N = n())
jam %>% subset(Month %in% c(3:6)) %>% group_by(Region, Month, Location, Sex.qPCR) %>% summarise(N = n())
jam %>% filter(!is.na(Month)) %>% filter(Region %in% "BEL") %>% nrow()  # 40
jam %>% filter(!is.na(Month), Month %nin% c(1:3,7,8,12)) %>% filter(Region %in% "BEL") %>% nrow()  # 33
jam %>% filter(Month %in% c(4:6,9:11)) %>% filter(Region %nin% "JAM") %>% nrow()  # 41
jam %>% filter(Month %in% c(4:6,9:11)) %>% group_by(Region, Month, Sex.qPCR) %>% summarise(N = n())
jam %>% filter(Month %in% c(1:3,12)) %>% group_by(Region, Month, Sex.qPCR) %>% summarise(N = n())

hbsc <- d[which(d$Membership.rad %in% "HBSC" | d$Membership.wgs %in% "HBSC"),]
hbsc %>% subset(Month %in% c(7,8)) %>% nrow()
hbsc %>% subset(Month %in% c(7,8)) %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)
hbsc %>% subset(Month %in% c(9:11)) %>% group_by(Region) %>% summarise(N = n()) %>% print(n = 50)
hbsc %>% group_by(Region, Sex.qPCR) %>% summarise(N = n()) %>% print(n = 50)

sle <- d[which(d$Membership.rad %in% "SLE" | d$Membership.wgs %in% "SLE"),]
sle %>% group_by(Month) %>% summarise(N = n())

reb <- d[which(d$Membership.wgs %in% "RES"),]
reb %>% group_by(Month) %>% summarise(N = n())

# By region
d %>% group_by(Region) %>% summarise(N = n())

stlawrence <- d[which(d$Region %in% "SLE"),]
stlawrence <- stlawrence %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
stlawrence %>% subset(Month %in% c(7,8)) %>% nrow()
stlawrence %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
stlawrence %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
stlawrence %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

cumbsound <- d[which(d$Region %in% "CS"),]
cumbsound <- cumbsound %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
cumbsound %>% subset(Month %in% c(7,8)) %>% nrow()
cumbsound %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
cumbsound %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
cumbsound %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

frobisher <- d[which(d$Region %in% "FRB"),]
frobisher <- frobisher %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
frobisher %>% subset(Month %in% c(7,8)) %>% nrow()
frobisher %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
frobisher %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
frobisher %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

ungava <- d[which(d$Region %in% "UNG"),]
ungava <- ungava %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
ungava %>% subset(Month %in% c(7,8)) %>% nrow()
ungava %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
ungava %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
ungava %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

Nhudstrait <- d[which(d$Region %in% "NHS"),]
Nhudstrait <- Nhudstrait %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
Nhudstrait %>% subset(Month %in% c(7,8)) %>% nrow()
Nhudstrait %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
Nhudstrait %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
Nhudstrait %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

Shudstrait <- d[which(d$Region %in% "SHS"),]
Shudstrait <- Shudstrait %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
Shudstrait %>% subset(Month %in% c(7,8)) %>% nrow()
Shudstrait %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
Shudstrait %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
Shudstrait %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

NEhudson <- d[which(d$Region %in% "NEH"),]
NEhudson <- NEhudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
NEhudson %>% subset(Month %in% c(7,8)) %>% nrow()
NEhudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
NEhudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
NEhudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

SEhudson <- d[which(d$Region %in% "SEH"),]
SEhudson <- SEhudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
SEhudson %>% subset(Month %in% c(7,8)) %>% nrow()
SEhudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
SEhudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
SEhudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

belcher <- d[which(d$Region %in% "BEL"),]
belcher <- belcher %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
belcher %>% subset(Month %in% c(7,8)) %>% nrow()
belcher %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
belcher %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
belcher %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

jamesbay <- d[which(d$Region %in% "JAM"),]
jamesbay <- jamesbay %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
jamesbay %>% subset(Month %in% c(7,8)) %>% nrow()
jamesbay %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
jamesbay %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
jamesbay %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

Whudson <- d[which(d$Region %in% c("SWH","NWH","NHB")),]
Whudson <- Whudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
Whudson %>% subset(Month %in% c(7,8)) %>% nrow()
Whudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
Whudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
Whudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

SWhudson <- d[which(d$Region %in% "SWH"),]
SWhudson <- SWhudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
SWhudson %>% subset(Month %in% c(7,8)) %>% nrow()
SWhudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
SWhudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
SWhudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

NWhudson <- d[which(d$Region %in% "NWH"),]
NWhudson <- NWhudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
NWhudson %>% subset(Month %in% c(7,8)) %>% nrow()
NWhudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
NWhudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
NWhudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

Nhudson <- d[which(d$Region %in% "NHB"),]
Nhudson <- Nhudson %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
Nhudson %>% subset(Month %in% c(7,8)) %>% nrow()
Nhudson %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
Nhudson %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
Nhudson %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())

resolute <- d[which(d$Region %in% "RES"),]
resolute <- resolute %>% mutate(Membership = ifelse(Marker %in% "ddRADseq", Membership.rad, Membership.wgs))
resolute %>% subset(Month %in% c(7,8)) %>% nrow()
resolute %>% subset(Month %in% c(7,8)) %>% group_by(Membership, Sex.qPCR) %>% summarise(N = n())
resolute %>% subset(Month %in% c(4:6)) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
resolute %>% subset(Month > 8) %>% group_by(Membership, Month, Sex.qPCR) %>% summarise(N = n())
