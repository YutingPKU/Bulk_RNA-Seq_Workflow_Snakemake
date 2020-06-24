#!/usr/bin/env Rscript
#-------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: May, 23, 2016
# @modified by Yuting Liu
# @modified date: Jun, 24, 2020
#--------------------



options(error = function() traceback(2))

library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(ggrepel)
suppressMessages(source('modules/scripts/supp_fns.R'))

#source("/mnt/cfce-stor1/home/mgc31/code/viperproject/modules/scripts/supp_fns.R")
#rpkmFile = "results/cufflinks/Cuff_Gene_Counts.filtered_top1000VarGenes.csv"
#metaFile = "data/metasheet.csv"
#pca_out_dir = "results/plots/"

pca_plot <- function(rpkmTable, annot, pca_out_dir) {
  rpkm.pca <- prcomp(t(rpkmTable), center = TRUE, scale. = TRUE)
  plot.var <- ggscreeplot(rpkm.pca)
  suppressMessages(ggsave(paste(pca_out_dir,"images/pca_plot_scree.png", sep=""), width = 1.4, height = 1.5, units = 'in'))
  

  res.pca <- rpkm.pca
  eig <- (res.pca$sdev)^2
  variance <- eig*100/sum(eig)
  cumvar <- cumsum(variance)
  eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                      cumvariance = cumvar)
  fontsize = 8
  linesize = 1
  
  all_plots <- list()
  for (ann in colnames(annot)){
    groups = as.character(annot[,ann])
    plotf <- function(gp){
      g <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = gp)) +
      geom_point(size =.8) + 
      #  geom_text(aes(label=colnames(rpkmTable)),hjust=.5, vjust=1.5, size = 1)+
      geom_label_repel(aes(label = colnames(rpkmTable)),
                                             box.padding   = 0.35, 
                                             point.padding = 0.5, size = 1)+
      xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
      ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
      ggtitle(paste0(" PCA log2 FPKM (", nrow(rpkmTable), "genes) "))+
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", size = linesize ))+
      theme(axis.text.x = element_text(size=fontsize),
            axis.text.y = element_text(size=fontsize),  
            axis.title.x = element_text(size=fontsize),
            axis.title.y = element_text(size=fontsize),
            axis.line = element_blank(), axis.ticks = element_line(size = .3),
            legend.position = "none", 
            #legend.text = element_text(size = fontsize),
            #legend.title  = element_text(size = fontsize),legend.key.size = unit(0.2, "cm"),
            plot.title = element_text(size=fontsize, hjust = 0.5))
      return(g)
    }
    g <- plotf(groups)
   
    all_plots[[ann]] = g
    suppressMessages(ggsave(paste(pca_out_dir, "images/pca_plot_", ann, ".png", sep=""), width = 4, height = 4, units = 'in', dpi = 1000))
  }

  pdf(paste(pca_out_dir, "pca_plot.pdf", sep=""), width = 4, height = 4)
  capture.output(print(c(all_plots,list(plot.var))))
  junk <- dev.off()
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile <- args[1]
metaFile <- args[2]
pca_out_dir <- args[3]

rpkmTable <- read.csv(rpkmFile, header=T, check.names=F,
                        row.names=1, stringsAsFactors=FALSE, dec='.')
annot <- read.csv(metaFile, sep=",", header=T, row.names=1,
                      stringsAsFactors=FALSE, check.names=F, comment.char='#')
if(any(grepl("comp_*", colnames(annot)))) {
  annot <- annot[, !grepl('comp_*', colnames(annot)), drop = F]
}
pca_plot(rpkmTable, annot, pca_out_dir)

