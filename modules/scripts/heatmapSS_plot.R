#!/usr/bin/env Rscript
#-------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: May, 23, 2016
# @modified by Yuting Liu
# @modified date: Jun, 24, 2020
#--------------------

## Load required packages
suppressMessages(library("gplots"))
suppressWarnings(suppressMessages(library("ComplexHeatmap")))
suppressMessages(library("circlize"))
suppressMessages(library("viridis"))
suppressMessages(library('dplyr'))
suppressMessages(source('modules/scripts/supp_fns.R'))

## Enable stack trace
#options(error = function() traceback(2))

heatmapSS_Spearman_plot <- function(rpkmTable,annot, ss_out_dir) {
    
    ## Read in and Log Transform Data
    Exp_data <-  rpkmTable
    #Exp_data <- log2(rpkmTable+1)
    #CHECK: DROP cols that are all 0s
    #LEN: on reverting, this line was causing the entire matrix to drop out
    #Exp_data <- Exp_data[, -(which(colSums(Exp_data) == 0))]

    ## Calc. spearman correlation
    cordata <- cor(Exp_data, method="spearman")

    # NOTES on clustering, not used for now
    # Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
    # Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
    rowdistance = dist(as.matrix(cordata), method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    
    ## make SS (sample-sample) heatmap
    ma_nolym <- max(cordata)
    mi_nolym <- min(cordata)
    my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)
    

    pdf(file = paste(ss_out_dir, "heatmapSS_Spearman_plot.pdf", sep=""), width = 2.5 + ncol(annot)*0.3, height = 1.8)

    ha1 <- make_complexHeatmap_annotation(annot)

    mapplot <-Heatmap(t(as.matrix(cordata)),
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #column_dend_height = unit(2, "cm"),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Sample Spearman Correlation",
                     column_title_gp = gpar(fontsize = 3.5),
                     #row_title = "Samples",
                     show_row_names = TRUE, show_column_names = FALSE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 3),
                     column_names_gp = gpar(fontsize = 0),
                     row_dend_width = unit(2.5,"mm"),
                     column_dend_height = unit(2.5, "mm"),
                     row_dend_gp = gpar(lwd = .35), column_dend_gp = gpar(lwd = .35),
                     cluster_rows = TRUE,
                     cluster_columns=TRUE,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param=list(title="corr", title_gp=gpar(fontsize=3.5), labels_gp=gpar(fontsize=3), legend_width = unit(.5,"mm")),
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
    draw(mapplot)
   # for(an in colnames(annot[1:ncol(annot)])) {
  #      decorate_annotation(an, {
  #          grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
  #          grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
  #      })
   # }
    dev.off()
    
    png(file=paste(ss_out_dir, "images/heatmapSS_Spearman_plot.png", sep=""), width = 2.5+ncol(annot)*0.3, height = 1.8, unit="in",res=1000)
    draw(mapplot)
   # for(an in colnames(annot[1:ncol(annot)])) {
  #      decorate_annotation(an, {
  #          grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
  #          grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
  #      })
  #  }
    dev.off()

    #WRITE output to file
    output<-as.matrix(cordata)
    output<-output[rowcluster$order, colcluster$order]
    write.table(output, file=paste(ss_out_dir, "heatmapSS_Spearman.txt",sep=""), quote=F, col.names = NA, sep="\t")
    
}
heatmapSS_Pearson_plot <- function(rpkmTable,annot, ss_out_dir) {
  
  ## Read in and Log Transform Data
  Exp_data <- rpkmTable
  #CHECK: DROP cols that are all 0s
  #LEN: on reverting, this line was causing the entire matrix to drop out
  #Exp_data <- Exp_data[, -(which(colSums(Exp_data) == 0))]
  
  ## Calc. spearman correlation
  cordata <- cor(Exp_data, method="pearson")
  
  # NOTES on clustering, not used for now
  # Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
  # Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
  rowdistance = dist(as.matrix(cordata), method = "euclidean")
  rowcluster = hclust(rowdistance, method = "ward.D2")
  coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
  colcluster = hclust(coldistance, method = "ward.D2")
  
  ## make SS (sample-sample) heatmap
  ma_nolym <- max(cordata)
  mi_nolym <- min(cordata)
  my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)
  
  pdf(file = paste(ss_out_dir, "heatmapSS_Pearson_plot.pdf", sep=""), width = 2.5+ncol(annot)*0.3, height = 1.8)
  
  ha1 <- make_complexHeatmap_annotation(annot)
  mapplot <-Heatmap(t(as.matrix(cordata)),
                    col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                    #column_dend_height = unit(2, "cm"),
                    #heatmap_legend_param = list(title = "exp. level"),
                    column_title = "Sample-Sample Pearson Correlation",
                    column_title_gp = gpar(fontsize = 3.5),
                    #row_title = "Samples",
                    show_row_names = TRUE, show_column_names = FALSE,
                    #row_names_max_width = unit(3, "mm"),
                    row_names_gp = gpar(fontsize = 3),
                    column_names_gp = gpar(fontsize = 0),
                    row_dend_width = unit(2.5,"mm"),
                    column_dend_height = unit(2.5, "mm"),
                    row_dend_gp = gpar(lwd = .35), column_dend_gp = gpar(lwd = .35),
                    cluster_rows = TRUE,
                    cluster_columns=TRUE,
                    show_heatmap_legend = TRUE,
                    heatmap_legend_param=list(title="corr", title_gp=gpar(fontsize=3.5), labels_gp=gpar(fontsize=3)),
                    
                    #row_dend_width = unit(5, "mm"),
                    #width=unit(60,"cm"),
                    top_annotation=ha1,
  )

  draw(mapplot)
  
  dev.off()
  
  png(file=paste(ss_out_dir, "images/heatmapSS_Pearson_plot.png", sep=""), width = 2.5+ncol(annot)*0.3, height = 1.8, unit="in",res=300)
  draw(mapplot)
  
  dev.off()
  
  #WRITE output to file
  output<-as.matrix(cordata)
  output<-output[rowcluster$order, colcluster$order]
  write.table(output, file=paste(ss_out_dir, "heatmapSS_Pearson.txt",sep=""), quote=F, col.names = NA, sep="\t")
  
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
ss_out_dir=args[3]

rpkmTable <- read.csv(rpkmFile, header=T, check.names=F, row.names=1, stringsAsFactors=FALSE, dec='.')

annot <- read.csv(annotFile, sep=",", header=T, row.names=1, stringsAsFactors=FALSE, check.names=F, comment.char='#')
if(any(grepl("comp_*", colnames(annot)))) {
  annot <- annot[, !grepl('Pair', colnames(annot)), drop = F]
  annot <- annot[, !grepl('comp_*', colnames(annot)), drop = F]
}

heatmapSS_Spearman_plot(rpkmTable,annot, ss_out_dir)
heatmapSS_Pearson_plot(rpkmTable,annot, ss_out_dir)
