#### violin plots ####
## create violinPlot function
violinPlot <- function(ctGenes, byFactor, factorOrder, groupLabel, 
                       extraLabel, dotSize, dotAlpha = 0.5){
  
  ## don't run if dataframe has less than four values
  if(nrow(ctGenes) > 4){
  
    #### set up plotting data ####
    ## account for kmeans.cluster column when in use
    if(byFactor != "kmeans.cluster"){
      idCols <- as.numeric(8)
      if("kmeans.cluster" %in% names(ctGenes)){
        ctGenes <- subset(ctGenes, select = -kmeans.cluster)
      }
    } else{
      idCols <- as.numeric(9)
      ## select specific subclusters if specified in factorOrder
      ctGenes <- subset(ctGenes, kmeans.cluster %in% factorOrder)
    }
    
    ## get p-values for differential expression for factors
    pvals <- NULL
    for(i in (idCols+1):ncol(ctGenes)){
      pvals <- c(pvals, kruskal.test(ctGenes[,i], factor(ctGenes[, byFactor]))$p.value)
    }
    pvals <- p.adjust(pvals, method = "BH")
    
    ctOrderGenes <- ctGenes[,c(1:idCols, (order(pvals)+idCols))]
    sigGeneNames <- names(ctGenes)[order(pvals)+idCols]
    sigGeneVals <- pvals[order(pvals)]
    sigGeneVals <- signif(sigGeneVals, digits=4)
    sigLength <- length(pvals[pvals < 0.05 & !is.na(pvals)])
    
    ## print differentially expressed genes
    print(paste("Differentially expressed genes between ", 
                groupLabel, " ", extraLabel, ":", sep=""), quote = F)
    print(paste(sigGeneNames[1:sigLength], sigGeneVals[1:sigLength], sep=": "), quote = F)
    
    ## gene set enrichment
    # gse(sigGeneNames, sigLength)
    
    ## melt for plotting
    ctMelt <- melt(ctOrderGenes, id.vars = 1:idCols, variable.name = "gene")
    names(ctMelt)[which(names(ctMelt)=="value")] <- "log2Ex"
    
    
    #### ggplot violins ####
#     genesToPlotGG1 <- sigGeneNames[1:6]
#     genesToPlotGG2 <- sigGeneNames[7:12]
#     genesToPlotGG3 <- sigGeneNames[13:18]
#     genesToPlotGG4 <- sigGeneNames[19:24]
    
    # genesToPlotGG <- list(genesToPlotGG1, genesToPlotGG2, genesToPlotGG3, genesToPlotGG4)
    
    genesToPlotGG <- list()
    ctMeltSubGG <- list()
    vplotList <- list()
    genesGG <- list()
    fillLabel <- list()
    for(i in 1:16){
      genesToPlotGG[[i]] <- sigGeneNames[(i*6 - 5):(i*6)]
      
      ctMeltSubGG[[i]] <- ctMelt[which(ctMelt$gene %in% genesToPlotGG[[i]]),]
      
      ctMeltSubGG[[i]] <- ctMeltSubGG[[i]][order(factor(ctMeltSubGG[[i]][, byFactor],
                                                        levels=factorOrder)),]
      
      genesGG[[i]] <- factor(ctMeltSubGG[[i]]$gene, levels=genesToPlotGG[[i]])
      fillLabel[[i]] <- factor(ctMeltSubGG[[i]][, byFactor], levels=factorOrder)
    }
    for(i in seq(1, 16, 4)){
#       vplot1 = ggplot(ctMeltSubGG[[i]], aes(genesGG[[i]], log2Ex, fill=fillLabel[[i]])) +
# #                            geom_rect(aes(xmin=genesGG[[i]],
# #                                           xmax=genesGG[[i]],
# #                                           ymin=-Inf,
# #                                           ymax=Inf,
# #                                           fill=fillLabel[[i]]))  +
#                            geom_violin(scale = "count", 
#                                        position=position_dodge(width = 0.75), 
#                                        bw = "nrd0",
#                                        adjust=3) +
#                            geom_point(size = dotSize, 
#                                       alpha = dotAlpha,
#                                       # shape = 21,    #colored circles
#                                       shape = 19, 
#                                       # stroke=1.35,   # stroke is for shape=1 circles
#                                       position=position_jitterdodge(jitter.width = 0.3, 
#                                                                     jitter.height = 0, 
#                                                                     dodge.width = 0.75)) +
#                            scale_fill_brewer(palette = "Set3") +
#                            # scale_color_brewer(palette = "Set3") +
#                            scale_y_continuous() +
#                            guides(fill=guide_legend(title=groupLabel)) +
#                            ggtitle(paste("most significant expression differences between ",
#                                          groupLabel, " ", extraLabel, " (plot #", i, ")", sep="")) +
#                            xlab("\ngenes") +
#                            ylab("log2(gene exp. / Gapdh exp.)\n") +
#                            theme_minimal() +
#                            theme(text=element_text(size=25),
#                                  panel.grid.minor=element_line(color="gray90"),
#                                  panel.grid.major=element_line(color="gray75", size=0.4),
#                                  panel.grid.major.x=element_blank(),
#                                  axis.title.x=element_text(vjust=-0.5),
#                                  plot.margin=unit(c(1,0,1,0),"cm"))
      vplot1 = ggplot(ctMeltSubGG[[i]], aes(genesGG[[i]], log2Ex, fill=fillLabel[[i]])) +
        geom_violin(scale = "width", 
                    position=position_dodge(width = 0.75), 
                    bw = "nrd0",
                    adjust=0.8) +
        geom_point(size = dotSize, 
                   alpha = dotAlpha, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title=groupLabel)) +
        ggtitle(paste("most significant expression differences between ",
                      groupLabel, " ", extraLabel, " (plot #", (i), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot2 = ggplot(ctMeltSubGG[[i+1]], aes(genesGG[[i+1]], log2Ex, fill=fillLabel[[i+1]])) +
                      geom_violin(scale = "width", 
                                  position=position_dodge(width = 0.75), 
                                  bw = "nrd0",
                                  adjust=0.8) +
                      geom_point(size = dotSize, 
                                 alpha = dotAlpha, 
                                 shape = 19, 
                                 # stroke=1.35,   # stroke is for shape=1 circles
                                 position=position_jitterdodge(jitter.width = 0.3, 
                                                               jitter.height = 0, 
                                                               dodge.width = 0.75)) +
                      scale_fill_brewer(palette = "Set3") +
                      scale_y_continuous() +
                      guides(fill=guide_legend(title=groupLabel)) +
                      ggtitle(paste("most significant expression differences between ",
                                    groupLabel, " ", extraLabel, " (plot #", (i+1), ")", sep="")) +
                      xlab("\ngenes") +
                      ylab("log2(gene exp. / Gapdh exp.)\n") +
                      theme_minimal() +
                      theme(text=element_text(size=25),
                            panel.grid.minor=element_line(color="gray90"),
                            panel.grid.major=element_line(color="gray75", size=0.4),
                            panel.grid.major.x=element_blank(),
                            axis.title.x=element_text(vjust=-0.5),
                            plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot3 = ggplot(ctMeltSubGG[[i+2]], aes(genesGG[[i+2]], log2Ex, fill=fillLabel[[i+2]])) +
                      geom_violin(scale = "width", 
                                  position=position_dodge(width = 0.75), 
                                  bw = "nrd0",
                                  adjust=0.8) +
                      geom_point(size = dotSize, 
                                 alpha = dotAlpha, 
                                 shape = 19, 
                                 # stroke=1.35,   # stroke is for shape=1 circles
                                 position=position_jitterdodge(jitter.width = 0.3, 
                                                               jitter.height = 0, 
                                                               dodge.width = 0.75)) +
                      scale_fill_brewer(palette = "Set3") +
                      scale_y_continuous() +
                      guides(fill=guide_legend(title=groupLabel)) +
                      ggtitle(paste("most significant expression differences between ",
                                    groupLabel, " ", extraLabel, " (plot #", (i+2), ")", sep="")) +
                      xlab("\ngenes") +
                      ylab("log2(gene exp. / Gapdh exp.)\n") +
                      theme_minimal() +
                      theme(text=element_text(size=25),
                            panel.grid.minor=element_line(color="gray90"),
                            panel.grid.major=element_line(color="gray75", size=0.4),
                            panel.grid.major.x=element_blank(),
                            axis.title.x=element_text(vjust=-0.5),
                            plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot4 = ggplot(ctMeltSubGG[[i+3]], aes(genesGG[[i+3]], log2Ex, fill=fillLabel[[i+3]])) +
                      geom_violin(scale = "width", 
                                  position=position_dodge(width = 0.75), 
                                  bw = "nrd0",
                                  adjust=0.8) +
                      geom_point(size = dotSize, 
                                 alpha = dotAlpha, 
                                 shape = 19, 
                                 # stroke=1.35,   # stroke is for shape=1 circles
                                 position=position_jitterdodge(jitter.width = 0.3, 
                                                               jitter.height = 0, 
                                                               dodge.width = 0.75)) +
                      scale_fill_brewer(palette = "Set3") +
                      scale_y_continuous() +
                      guides(fill=guide_legend(title=groupLabel)) +
                      ggtitle(paste("most significant expression differences between ",
                                    groupLabel, " ", extraLabel, " (plot #", (i+3), ")", sep="")) +
                      xlab("\ngenes") +
                      ylab("log2(gene exp. / Gapdh exp.)\n") +
                      theme_minimal() +
                      theme(text=element_text(size=25),
                            panel.grid.minor=element_line(color="gray90"),
                            panel.grid.major=element_line(color="gray75", size=0.4),
                            panel.grid.major.x=element_blank(),
                            axis.title.x=element_text(vjust=-0.5),
                            plot.margin=unit(c(1,0,1,0),"cm"))
      
      grid.arrange(vplot1, vplot2, vplot3, vplot4, ncol=1, nrow=4)
    } ## ggplot violins

    
    #### trellis violins ####
#     genesToPlotTrellis <- names(ctOrderGenes[, (idCols+1):ncol(ctOrderGenes)])
#     
#     ctMeltSubTrellis <- ctMelt[which(ctMelt$gene %in% genesToPlotTrellis),]
#     
#     ctMeltSubTrellis <- ctMeltSubTrellis[order(factor(ctMeltSubTrellis[, byFactor], 
#                                                       levels=factorOrder)),]
#     
#     genesTrellis <- factor(ctMeltSubTrellis$gene, 
#                             levels=names(ctOrderGenes[, (idCols+1):ncol(ctOrderGenes)]))
#     
#     windowNumber <- length(levels(genesTrellis))
#     windowVert <- ceiling(windowNumber/4)
#     
#     lattice.options(as.table=TRUE)
#     dots <- stripplot(log2Ex ~ ctMeltSubTrellis[, byFactor] | genesTrellis, 
#                       data = ctMeltSubTrellis,
#                       jitter.data = TRUE, 
#                       alpha = 0.6,
#                       main = paste("Expression ordered by largest difference between ", groupLabel, " ", extraLabel, sep=""),
#                       xlab = groupLabel,
#                       ylab = "log2(gene exp. / Gapdh exp.)", 
#                       horizontal=F,
#                       as.table=TRUE,
#                       pch=1,
#                       scales=list(x=list(alternating=TRUE, cex=1.1), y=list(alternating=TRUE)),
#                       layout=(c(4,windowVert)))
#     
#     violin <- bwplot(log2Ex ~ ctMeltSubTrellis[, byFactor] | genesTrellis,
#                      data=ctMeltSubTrellis,
#                      outer=TRUE, 
#                      as.table=TRUE, 
#                      horizontal=FALSE,
#                      col='transparent',
#                      panel=panel.violin,
#                      layout=(c(4,windowVert)))
#     
#     dots + as.layer(violin)
  }
}