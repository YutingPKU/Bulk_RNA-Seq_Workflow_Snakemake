## Load required packages
suppressMessages(library("gplots"))
suppressWarnings(suppressMessages(library("ComplexHeatmap")))
suppressMessages(library("circlize"))
suppressMessages(library("viridis"))
suppressMessages(library("dplyr"))
suppressMessages(source('modules/scripts/supp_fns.R'))

## Enable stack trace
#options(error = function() traceback(2))

heatmapSF_Spearman_plot <- function(rpkmTable,annot, num_kmeans_clust, sf_out_dir) {
    
    ## Read in and Log Transform Data
    Exp_data = log2(rpkmTable+1)
    #CHECK: DROP cols that are all 0s
    #LEN: on reverting, this line was causing the entire matrix to drop out
    #Exp_data <- Exp_data[, -(which(colSums(Exp_data) == 0))]
    
    ## Calc. spearman correlation and use values for column clustering before any other alterations
    cordata <- cor(Exp_data, method="spearman")

    coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    
    ## Make SF (sample-feature) heatmap
    Exp_data <- apply(Exp_data,1,function(x) zscore(x))

    ## Set breaks for data for color scale
    ma_nolym <- max(Exp_data)
    mi_nolym <- min(Exp_data)
    my.breaks_nolym<-c(-3,seq(-2.5,2.5,length.out=99),3)

    ## Data needs to be transposed for heatmap
    Exp_data = t(as.matrix(Exp_data))

    ## Make annotation bars
    ha1 <- make_complexHeatmap_annotation(annot)

    ## Turn on rownames if less than 100 genes
    row_name_param = FALSE
    if (nrow(Exp_data) <= 100) {row_name_param = TRUE}

    ## Determine type of plot, and plot
    if (is.numeric(num_kmeans_clust) == TRUE) {kmparam = num_kmeans_clust}
    if (is.character(num_kmeans_clust) == TRUE) {kmparam = as.numeric(unlist(strsplit(num_kmeans_clust,",")))}
    
    pdf(file = paste(sf_out_dir, "heatmapSF_Spearman_plot.pdf", sep=""), width=11,height=8.5) 
    
    png_count = 0
    for (i in 1:length(kmparam)) {

        ## If kmparam is 0, use hierarichical
        if (kmparam[i] == 0) {
            rowdistance = dist(Exp_data, method = "euclidean")
            rowcluster = hclust(rowdistance, method = "ward.D2")
            rowclusterparam = rowcluster
            hmdata = Exp_data
            column_title_param = "Sample-Feature Hierarchical Clustering"
        }

        ## If kmparam is not 0, use kmeans
        if (kmparam[i] != 0) {
            #if (kmparam[i] == 1) {kmparam[i] = 2}
            km1 = kmeans(Exp_data, centers=kmparam[i])
            kmclust = km1$cluster
            kmclustsort = sort(kmclust)
            ind = match(names(kmclustsort), rownames(Exp_data))
            hmdata = Exp_data[ind,]
            rowclusterparam = FALSE
            column_title_param = paste("Sample-Feature Kmeans ", kmparam[i], " Clustering", sep="")
        }        
        
        mapplot <- Heatmap(hmdata,
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = column_title_param,
                     show_row_names = row_name_param, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 8),
                     cluster_rows = rowclusterparam,
                     cluster_columns = colcluster,
                     show_heatmap_legend = FALSE,
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
        
        ## First drawing into png
        png_count = png_count+1
        png(file=paste(sf_out_dir, "images/heatmapSF_Spearman_",png_count,"_plot.png",sep=""), width = 8, height = 8, unit="in",res=300)
        draw(mapplot)
        for(an in colnames(annot[1:ncol(annot)])) {
            decorate_annotation(an,
              {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
              grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
              })
        }
        junk <- dev.off()

        ## Repeated to get into the pdf
        draw(mapplot)
        for(an in colnames(annot[1:ncol(annot)])) {
            decorate_annotation(an,
              {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
              grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
              })
        }
        if (i == 1) {
            if (kmparam[1] == 0) {
                output<-Exp_data
                output<-output[unlist(row_order(mapplot)), unlist(column_order(mapplot))]
                write.table(output, file=paste(sf_out_dir, "heatmapSF_Spearman.txt",sep=""), quote=F, col.names = NA, sep="\t")
            }
            if (kmparam[1] != 0) {
                output = cbind(hmdata,kmclustsort)
                output = output[,unlist(column_order(mapplot))]
                write.table(output, file=paste(sf_out_dir, "heatmapSF_Spearman.txt",sep=""), quote=F, col.names = NA, sep="\t")
            }
        }
    }
    junk <- dev.off()        
        
}

heatmapSF_Pearson_plot <- function(rpkmTable,annot, num_kmeans_clust, sf_out_dir) {
  
  ## Read in and Log Transform Data
  Exp_data = log2(rpkmTable+1)
  #CHECK: DROP cols that are all 0s
  #LEN: on reverting, this line was causing the entire matrix to drop out
  #Exp_data <- Exp_data[, -(which(colSums(Exp_data) == 0))]
  
  ## Calc. spearman correlation and use values for column clustering before any other alterations
  cordata <- cor(Exp_data, method="pearson")
  
  coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
  colcluster = hclust(coldistance, method = "ward.D2")
  
  ## Make SF (sample-feature) heatmap
  Exp_data <- apply(Exp_data,1,function(x) zscore(x))
  
  ## Set breaks for data for color scale
  ma_nolym <- max(Exp_data)
  mi_nolym <- min(Exp_data)
  my.breaks_nolym<-c(-3,seq(-2.5,2.5,length.out=99),3)
  
  ## Data needs to be transposed for heatmap
  Exp_data = t(as.matrix(Exp_data))
  
  ## Make annotation bars
  ha1 <- make_complexHeatmap_annotation(annot)
  
  ## Turn on rownames if less than 100 genes
  row_name_param = FALSE
  if (nrow(Exp_data) <= 100) {row_name_param = TRUE}
  
  ## Determine type of plot, and plot
  if (is.numeric(num_kmeans_clust) == TRUE) {kmparam = num_kmeans_clust}
  if (is.character(num_kmeans_clust) == TRUE) {kmparam = as.numeric(unlist(strsplit(num_kmeans_clust,",")))}
  
  pdf(file = paste(sf_out_dir, "heatmapSF_Pearson_plot.pdf", sep=""), width=11,height=8.5) 
  
  png_count = 0
  for (i in 1:length(kmparam)) {
    
    ## If kmparam is 0, use hierarichical
    if (kmparam[i] == 0) {
      rowdistance = dist(Exp_data, method = "euclidean")
      rowcluster = hclust(rowdistance, method = "ward.D2")
      rowclusterparam = rowcluster
      hmdata = Exp_data
      column_title_param = "Sample-Feature Hierarchical Clustering"
    }
    
    ## If kmparam is not 0, use kmeans
    if (kmparam[i] != 0) {
      #if (kmparam[i] == 1) {kmparam[i] = 2}
      km1 = kmeans(Exp_data, centers=kmparam[i])
      kmclust = km1$cluster
      kmclustsort = sort(kmclust)
      ind = match(names(kmclustsort), rownames(Exp_data))
      hmdata = Exp_data[ind,]
      rowclusterparam = FALSE
      column_title_param = paste("Sample-Feature Kmeans ", kmparam[i], " Clustering", sep="")
    }        
    
    mapplot <- Heatmap(hmdata,
                       col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                       #heatmap_legend_param = list(title = "exp. level"),
                       column_title = column_title_param,
                       show_row_names = row_name_param, show_column_names = TRUE,
                       #row_names_max_width = unit(3, "mm"),
                       row_names_gp = gpar(fontsize = 6),
                       column_names_gp = gpar(fontsize = 8),
                       cluster_rows = rowclusterparam,
                       cluster_columns = colcluster,
                       show_heatmap_legend = FALSE,
                       #row_dend_width = unit(5, "mm"),
                       #width=unit(60,"cm"),
                       top_annotation=ha1,
    )
    
    ## First drawing into png
    png_count = png_count+1
    png(file=paste(sf_out_dir, "images/heatmapSF_Pearson",png_count,"_plot.png",sep=""), width = 8, height = 8, unit="in",res=300)
    draw(mapplot)
    for(an in colnames(annot[1:ncol(annot)])) {
      decorate_annotation(an,
                          {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
                            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
                          })
    }
    junk <- dev.off()
    
    ## Repeated to get into the pdf
    draw(mapplot)
    for(an in colnames(annot[1:ncol(annot)])) {
      decorate_annotation(an,
                          {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
                            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
                          })
    }
    if (i == 1) {
      if (kmparam[1] == 0) {
        output<-Exp_data
        output<-output[unlist(row_order(mapplot)), unlist(column_order(mapplot))]
        write.table(output, file=paste(sf_out_dir, "heatmapSF_Pearson.txt",sep=""), quote=F, col.names = NA, sep="\t")
      }
      if (kmparam[1] != 0) {
        output = cbind(hmdata,kmclustsort)
        output = output[,unlist(column_order(mapplot))]
        write.table(output, file=paste(sf_out_dir, "heatmapSF_Pearson.txt",sep=""), quote=F, col.names = NA, sep="\t")
      }
    }
  }
  junk <- dev.off()        
  
}


## Read in arguments
args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
num_kmeans_clust=args[3]
sf_out_dir=args[4]

rpkmTable <- read.csv(rpkmFile, header=T, check.names=F, row.names=1, stringsAsFactors=FALSE, dec='.')

annot <- read.csv(annotFile, sep=",", header=T, row.names=1, stringsAsFactors=FALSE, check.names=F, comment.char='#')
if(any(grepl("comp_*", colnames(annot)))) {
  annot <- annot[, !grepl('Pair', colnames(annot)), drop = F]
  annot <- annot[, !grepl('comp_*', colnames(annot)), drop = F]
}
## Run the function
heatmapSF_Spearman_plot(rpkmTable,annot, num_kmeans_clust, sf_out_dir)
heatmapSF_Pearson_plot(rpkmTable,annot, num_kmeans_clust, sf_out_dir)