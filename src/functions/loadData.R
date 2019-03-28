#### Load Data ####
loadData <- function(){
  
  # fileDir <- "~/Documents/teyton_handover/single_cell_macrophage/data/"
  fileDir <- paste(baseDir, "data/", sep="")
  
  fileList <- list.files(fileDir)[grep("meta.csv", list.files(fileDir), invert = T)]
  
  dataFiles <- paste(fileDir, fileList, sep="/")
  
  ## check for metadata files
  # fileDir <- dirname(dataFiles[1])
  metaFiles <- list.files(fileDir)[grep(".meta.csv", list.files(fileDir))]
  
  # if(length(metaFiles) == 0){
  #   showModal(modalDialog(
  #     title = "No Metadata Files Found",
  #     "Please create metadata file(s)",
  #     easyClose = TRUE
  #   ))
  # }
  
  if(length(metaFiles) > 0){
    ## read in data files and save data frames in dataList
    dataList <- list()
    metaDataList <- list()
    for(file in dataFiles){
      ## get file name
      fileName <- gsub(".csv", "", gsub(".*/", "", file))
      ## read data and save
      dataList[[fileName]] <- as.data.frame(read.csv(file, 
                                                     header = TRUE, 
                                                     stringsAsFactors = FALSE, 
                                                     skip = 11), 
                                            stringsAsFactors = F)
      
      ## read metadata
      if(length(metaFiles) == length(dataFiles)){
        metaDataFile <- gsub(".csv", ".meta.csv", file)
      }
      
      if(length(metaFiles) == 1){
        metaDataFile <- file.path(fileDir, metaFiles)
      }
      
      if(file.exists(metaDataFile)){
        metaDataList[[fileName]] <- as.data.frame(read.csv(metaDataFile, 
                                                           header = TRUE, 
                                                           stringsAsFactors = FALSE), 
                                                  stringsAsFactors = F)
        # change all metadata columns to character vectors (for now)
        metaDataList[[fileName]][] <- lapply(metaDataList[[fileName]], as.character)
      }
    }
  }
  
  ## format data frames
  ctTableList <- list()
  for(plate in names(dataList)){
    ### format data frames
    ## only keep relevant columns
    dataList[[plate]] <- dataList[[plate]][, c("Name", "Name.1", "Value", "Call")]
    ## change column names
    names(dataList[[plate]]) <- c("cellID", "gene", "ct", "Call")
    
    ### merge metadata
    if(length(metaDataList) > 0){
      ctTableList[[plate]] <- merge(dataList[[plate]], 
                                    metaDataList[[plate]], 
                                    by= "cellID", 
                                    all.x = T)
    }
    
    ## check for missing data and metadata
    
    ## create plate column
    ctTableList[[plate]]$plateID <- plate
    
    ### calculate log2Expression
    ## basic calculations based on recommended LoD from singular analysis toolkit
    ctTableList[[plate]]$log2Ex <- (24 - ctTableList[[plate]]$ct)
    ## set negative values to zero
    ctTableList[[plate]][(ctTableList[[plate]]$log2Ex < 0), "log2Ex"] <- as.numeric(0)
    ## set failed runs to Zero
    ctTableList[[plate]][(
      ctTableList[[plate]]$Call=="Fail" & ctTableList[[plate]]$ct!=999
    ), "log2Ex"] <- as.numeric(0)
  }
  
  ## concatenate data frames
  # ctTable <- do.call("rbind", ctTableList)
  ctTable <- rbind.fill(ctTableList)
  
  ### create final format
  ## remove "ct" and "Call" columns
  ctTable <- subset(ctTable, select = -c(ct, Call))
  
  ## rename duplicate cellIDs
  # just add plateID to cellIDs
  # assuming that the biomark software makes unique cellIDs for each plate
  ctTable$cellID <- paste(ctTable$cellID, ctTable$plateID, sep = "-")
  
  ## melt/cast to reformat
  ctCast <- recast(ctTable, ... ~ gene, measure.var = "log2Ex")
  
  ctCast <- subset(ctCast, select = -variable)
  
  ## Do not require complete cases (allows NAs)
  #ctCast <- ctCast[complete.cases(ctCast), ]
  
  
  return(ctCast)
}