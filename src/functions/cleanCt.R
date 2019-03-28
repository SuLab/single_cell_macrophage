#### clean Ct table ####
require(reshape2)

cleanCt <- function(ctRaw, summaryOutput = FALSE, cumExpCutoff = 0, cumHist = F, normGene = NULL){
  
  ## set gene for normalization
  if(is.null(normGene)){
    normGene <- "Gapdh"
  }
  
  #### calculate log2Expression ####
  ## basic calculations based on recommended LoD from singular analysis toolkit
  ctRaw$log2Ex <- (24 - ctRaw$ct)
  ## set negative values to zero
  ctRaw[(ctRaw$log2Ex < 0), "log2Ex"] <- as.numeric(0)
  ## set failed runs to Zero
  ctRaw[(ctRaw$Call=="Fail" & ctRaw$ct!=999), "log2Ex"] <- as.numeric(0)
  
  
  #### reformat data ####
  ## keep relevant columns
  ctByGene <- ctRaw[,c("mouse", "cellType", "gene", "log2Ex")]
  
  ## melt/cast to reformat
  ctCast <- recast(ctByGene, mouse + cellType ~ gene)
  
  ## change NAs to 0
  ctCast[is.na(ctCast)] <- as.numeric(0)
  
  ## remove 'blank' genes
  drops <- c("blank1", "blank2", "blank3")
  ctCast <- ctCast[ , !(names(ctCast) %in% drops)]
  
  #### remove outliers and normalize ####
  ## remove cells with no expression of any gene
  cellsTotal <- nrow(ctCast)
  cellsWithNoExp <- nrow(ctCast[which(rowSums(ctCast[,9:ncol(ctCast)], na.rm=T) == 0),])
  
  ctCast <- ctCast[which(rowSums(ctCast[,9:ncol(ctCast)], na.rm=T) > 0),]
  
  if(summaryOutput==T){
    cat(paste("No expression detected in ", cellsWithNoExp, "/", cellsTotal, " cells", sep=""))
  }
  
  if(normGene != "none"){
    ## remove cells with low gapdh expression
    if(summaryOutput==T){
      par(mar=c(50, 50, 50, 50))
      hist <- ggplot(ctCast, aes(ctCast[, normGene])) + 
        geom_histogram(colour = "blue", fill = "cyan", binwidth = 1) + 
        scale_x_continuous(breaks = seq(0, 16, 1), minor_breaks = NULL) +
        ggtitle(paste("histogram of raw ", normGene, " expression values per cell", sep="")) +
        xlab(paste(normGene, " log2Ex", sep="")) +
        theme(text=element_text(size=25),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(15,2,15,2),"cm"))
      print(hist)
      par(mar=c(5, 5, 5, 5))
    }
    
    totalCells <- nrow(ctCast)
    # ctCast <- subset(ctCast, normGene > 0)
    
    medianNormGene <- median(ctCast[, normGene])
    ctCast <- ctCast[ctCast[,normGene] > 8, ]
    # ctCast <- ctCast[ctCast[,normGene] > (medianNormGene/4), ]
    # ctCast <- subset(ctCast, normGene > (mediannormGene/2))
    removedCells <- totalCells - nrow(ctCast)

    if(summaryOutput==T){
      cat(paste(removedCells, "/", totalCells, " cells were removed due to low ", normGene, " expression (< ",
                signif((medianNormGene/4), digits=3), ")", sep=""))
    }
    
    #### MANUAL TESTING ####
    # medianNormGene <- 10
    # ctCast <- ctCast[ctCast[,normGene] > (medianNormGene), ]
    # # ctCast <- subset(ctCast, normGene > (mediannormGene/2))
    # removedCells <- totalCells - nrow(ctCast)
    # 
    # if(summaryOutput==T){
    #   cat(paste(removedCells, "/", totalCells, " cells were removed due to low ", normGene, " expression (< ",
    #             signif((medianNormGene), digits=3), ")", sep=""))
    # }
    
    #### Normalize by normGene Expression ####
    ## TESTING normalize normGene to 10 first
  #   gapdhNorm <- 10 - median(ctCast$normGene)
  #   ctCast$normGene <- ctCast$normGene + gapdhNorm
    ## assuming consistent normGene expression
    ctNorm <- ctCast[,9:ncol(ctCast)] - ctCast[, normGene]
    ## scale back to 0 expression reference
    # ctNorm <- ctNorm + median(ctCast$normGene)
    ## NO normGene Normalization ##
    # ctNorm <- ctCast[,9:ncol(ctCast)]
  
    #### Normalize zero expression to zero ####
  #   ctCastRawVals <- ctCast[,9:ncol(ctCast)]
  #   
  #   zeroMeanList <- list()
  #   for(gene in names(ctCastRawVals)){
  #     zeroMeanList[gene] <- mean(ctCastRawVals[which(ctCastRawVals[, gene] == 0), gene] - ctNorm[which(ctCastRawVals[, gene] == 0), gene])
  #   }
  #   
  #   zeroMeanList <- lapply(zeroMeanList, function(x) ifelse(is.nan(x),as.numeric(0),x))
  #   
  #   for(gene in names(ctNorm)){
  #     ctNorm[, gene] <- ctNorm[, gene] + zeroMeanList[[gene]]
  #   }
      
    ## REMOVED FOR TESTING
    ## center data at median
    for(i in 1:ncol(ctNorm)){
      ctNorm[, i] <- ctNorm[, i] - median(ctNorm[, i])
    }
    
    ctNorm <- cbind(ctCast[,1:8], ctNorm)
  }
  
  if(normGene == "none"){
    ctNorm <- ctCast
  }
  
  ## calculate total expression and remove outliers
  ctNorm$ctSum <- rowSums(ctNorm[,9:ncol(ctNorm)])
  
  if(summaryOutput==T){
    par(mar=c(50, 50, 50, 50))
    
    if(cumHist == T){
      hist <- ggplot(ctNorm, aes(ctSum)) + 
        geom_histogram(colour = "blue", fill = "cyan", binwidth = 5) + 
        # scale_x_continuous(breaks = seq(-300, 350, 50), minor_breaks = seq(-290, 340, 10)) +
        ggtitle("histogram of cumulative expression values per cell") +
        xlab("sum of normalized log2Ex values") +
        theme(text=element_text(size=25),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(15,2,15,2),"cm"))
      print(hist)
    }
  }
  
  if(cumExpCutoff!=F){
    totalCells <- nrow(ctNorm)
    ctNorm <- subset(ctNorm, ctSum > cumExpCutoff)
    remCells <- totalCells - nrow(ctNorm)
      
    if(summaryOutput==T){
      cat(paste(remCells, " / ", totalCells, " cells were removed due to low gene expression (< ", cumExpCutoff, ")", sep=""))
      
      if(cumHist == T){
        hist <- ggplot(ctNorm, aes(ctSum)) + 
          geom_histogram(colour = "blue", fill = "cyan", binwidth = 5) + 
          # scale_x_continuous(breaks = seq(cumExpCutoff, 350, 50), minor_breaks = seq(cumExpCutoff, 340, 10)) +
          ggtitle("histogram of cumulative expression values per cell after outlier removal") +
          xlab("sum of normalized log2Ex values") +
          theme(text=element_text(size=25),
                axis.title.x=element_text(vjust=-0.5),
                plot.margin=unit(c(15,2,15,2),"cm"))
        print(hist)
      }
    }
  }
  
  ## delete ctSum column
  ctNorm$ctSum <- NULL
  
  par(mar=c(5, 5, 5, 5))
  
  return(ctNorm)
}