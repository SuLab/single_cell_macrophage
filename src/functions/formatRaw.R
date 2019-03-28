#### format raw data ####
formatRaw <- function(){
  
  ## directory names
  inDir <- paste(baseDir, "data/raw/", sep="")
  outDir <- paste(baseDir, "data/formatted/", sep="")
  
  ## get file names
  filesObscurin <- list.files(paste(inDir, "obscurin", sep=""))
  filesObscurin <- paste("obscurin/", filesObscurin, sep="")
  
  files9D9Q <- list.files(paste(inDir, "9D_9Q", sep=""))
  files9D9Q <- paste("9D_9Q/", files9D9Q, sep="")
  
  filesInsulin <- list.files(paste(inDir, "insulin", sep=""))
  filesInsulin <- paste("insulin/", filesInsulin, sep="")
  
  filesPd1Icos <- list.files(paste(inDir, "pd1_icos", sep=""))
  filesPd1Icos <- paste("pd1_icos/", filesPd1Icos, sep="")
  
  filesGfp <- list.files(paste(inDir, "gfp", sep=""))
  filesGfp <- paste("gfp/", filesGfp, sep="")
  
  filesHuman <- list.files(paste(inDir, "human", sep=""))
  filesHuman <- paste("human/", filesHuman, sep="")
  
  filesCholx <- list.files(paste(inDir, "cholx", sep=""))
  filesCholx <- paste("cholx/", filesCholx, sep="")
  
  filesMacrophage <- list.files(paste(inDir, "macrophage", sep=""))
  filesMacrophage <- paste("macrophage/", filesMacrophage, sep="")
  
  filesB57 <- list.files(paste(inDir, "b57", sep=""))
  filesB57 <- paste("b57/", filesB57, sep="")
  
  ## concatenate file names
  files <- c(filesObscurin, files9D9Q, filesInsulin, filesPd1Icos, filesGfp, filesHuman, filesCholx, filesMacrophage, filesB57)

  ## run through files and format if formatted (.trimmed) file does not exist
  for(file in files){
    if(file.exists(paste(outDir, file, ".trimmed", sep="")) == F){
      print(file)
      ## convert to unix encoding
      system(paste("dos2unix ", inDir, file, sep=""))
      
      ## remove header
      system(paste("sed '1,11d' ", inDir, file, " > ", outDir, file, ".trimmed", sep=""))
    }
  }
}