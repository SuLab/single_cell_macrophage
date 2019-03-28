#### load macrophage data ####
dataLoadMacrophage <- function(){
  
  ## set data directory
  dataDir <- paste(baseDir, "data/formatted/macrophage/", sep="")
  
  ### read in data
  ## age: 7 weeks
  ## tissue: islets
  ## probes: macrophage
  ## sort: 1
  ## mice: 3 mice pooled (called mouse 1)
  ctTable1 <- read.csv(paste(dataDir, "macrophage_1_1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTable2 <- read.csv(paste(dataDir, "macrophage_1_2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTableSort1 <- rbind (ctTable1, ctTable2)
  
  ctTableSort1$mouse <- "1"
  
  
  ## sort: 3
  ## mice: 4, 5, and 6
  ctTable3 <- read.csv(paste(dataDir, "macrophage_3_1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTable4 <- read.csv(paste(dataDir, "macrophage_3_2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTableSort2 <- rbind (ctTable3, ctTable4)
  
  ctTableSort2$mouse <- substr(ctTableSort2$Name, 4, 4)
  
  
  ## sort: 4
  ## mice: 7, 8, and 9
  
  ctTable5 <- read.csv(paste(dataDir, "macrophage_4_1.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTable6 <- read.csv(paste(dataDir, "macrophage_4_2.csv.trimmed", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE)
  
  ctTableSort3 <- rbind (ctTable5, ctTable6)
  
  ctTableSort3$mouse <- substr(ctTableSort3$Name, 4, 4)
  
  ctTable <- rbind(ctTableSort1, ctTableSort2, ctTableSort3)
  
  ctTable <- ctTable[grep("mac", ctTable$Name), ]
  
  ## change column names
  names(ctTable)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## create other columns
  # ctTable$cellSource <- "I"
  # 
  # ctTable$probe <- "macrophage"
  # 
  # ctTable$age <- "7"
  
  return(ctTable)
  
}