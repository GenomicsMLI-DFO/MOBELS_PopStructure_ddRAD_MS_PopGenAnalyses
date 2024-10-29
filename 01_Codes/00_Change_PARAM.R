# Info --------------------------------------------------------------------
# 
# Authors: Audrey Bourret, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Institut Maurice Lamontagne 
# Date: 2021-01-12
# 
# Overview: This file allowed to set/modify folder structure
# 
#


# Library -----------------------------------------------------------------

# Internal functions
if(length(list.files("./01_Codes/Functions"))>=1){ 
  for(i in 1:length( list.files("./01_Codes/Functions") )){
    source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
  }
}

PARAM <- data.frame(x = character(), value = character())


# Specific functions ------------------------------------------------------

"add.param<-" <- function(PARAM, ..., value){
  PARAM <- rbind(PARAM, value)
  names(PARAM) <- c("x", "value")
  PARAM$x     <- as.character(PARAM$x)
  PARAM$value <- as.character(PARAM$value)
  PARAM
  
}


# Path --------------------------------------------------------------------

# 00_Data
add.param(PARAM) <- c("info.path", "../MOBELS_PopStructure_ddRAD_MS_PopGenAnalyses/00_Data/00_FileInfos")
add.param(PARAM) <- c("stacks.ref.path.bringloe.all", "../MOBELS_PopStructure_ddRAD_MS_PopGenAnalyses/00_Data/05b_Stacks.ref/01_Bringloe_22xii22/00_all") # External disk
add.param(PARAM) <- c("stacks.ref.path.bringloe.popgen", "../MOBELS_PopStructure_ddRAD_MS_PopGenAnalyses/00_Data/05b_Stacks.ref/01_Bringloe_22xii22/01_popgen") # External disk
add.param(PARAM) <- c("filter.ref.path.bringloe", "./00_Data/01_Filtering.ref/01_Bringloe")


# Results - Stacks
add.param(PARAM) <- c("stacks.ref.log.bringloe.all", "../MOBELS_PopStructure_ddRAD_MS_PopGenAnalyses/02_Results/00_Stacks/05b_Stacks.ref/01_Bringloe_22xii22/00_all")
add.param(PARAM) <- c("stacks.ref.log.bringloe.popgen", "../MOBELS_PopStructure_ddRAD_MS_PopGenAnalyses/02_Results/00_Stacks/05b_Stacks.ref/01_Bringloe_22xii22/01_popgen")


# Update file before continue to ensure next section will do OK!
if(dir.exists("./00_Data/00_FileInfos") == F){
  dir.create("./00_Data/00_FileInfos")
}

write.csv2(PARAM, file = file.path("./00_Data/00_FileInfos", "Options.csv"), row.names=F)




