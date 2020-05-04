#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 25, 2016
#---------------------------------

#suppressMessages(library(argparse))
if( is.element("Seurat", installed.packages())){
    suppressMessages(library(Seurat))
} else {
    suppressMessages(require("devtools"))
    source('http://bioconductor.org/biocLite.R')
    install_github("Storeylab/lfa")
    install_github("satijalab/seurat")
    suppressMessages(require(Seurat))
}

options(error = function() traceback(2))

runSeurat <- function( matrix_file, metasheet, out_dir ) {
    sc.data <- read.csv( matrix_file, header=TRUE, row.names=1 )
    sc.data <- log( sc.data + 1 )
    sc_obj <- new( "seurat", raw.data=sc.data )

    #------------------------------------
    # loop through annotations
    #------------------------------------
    ann <- read.csv(metasheet, header=T, row.names=1, as.is=T, comment.char='#')
    ann <- ann[,!grepl("comp_*", colnames( ann )), drop=F]

    for( index in range(1,length(colnames(ann))) ) {
        write(paste( "Processing ", colnames(ann)[index], sep=""), stderr() )
        dir.create( file.path(out_dir,colnames(ann)[index]), showWarnings=F)
        sc <- Setup( sc_obj, project="Viper Single Cell Analysis", min.cells = 3, names.field = index, names.delim = "_", min.genes = 1000, is.expr=1, )
        png( paste(out_dir, "/", colnames(ann)[index], "/mean_var_plot.png", sep=""), width = 8, height = 8, unit="in",res=300 )
        sc <- MeanVarPlot( sc, y.cutoff = 2, x.low.cutoff = 1, fxn.x = expMean, fxn.y = logVarDivMean )
        dev.off()
        png( paste(out_dir, "/", colnames(ann)[index], "/linear_PCA.png", sep=""), width = 8, height = 8, unit="in",res=300)
        sc <- PCA( sc, do.print=FALSE )
        PCAPlot( sc, 1, 2, pt.size = 2 )
        dev.off()
    
        png( paste(out_dir, "/", colnames(ann)[index], "/non_linear_tSNE_PCA.png", sep=""), width = 8, height = 8, unit="in",res=300 )
        sc <- RunTSNE( sc, dims.use = 1:length(unique(ann[,index])), max_iter=2000 )
        TSNEPlot( sc )
        dev.off()

        png( paste(out_dir, "/", colnames(ann)[index], "/viz_PCA.png", sep=""), width = 8, height = 8, unit="in",res=300 )
        VizPCA(sc,1:length(unique(ann[,index])))
        dev.off()
    } 
}

parse_args <- function() {
    parser <- argparse::ArgumentParser(description="Takes TPM matrix file and runs seurat tSNE methods")
    parser$add_argument('-m', '--matrix_file', help="TPM matrix file", type="character", nargs=1)
    args <- parser$parse_args()
    return (args)
}

#args <- parse_args()
args <- commandArgs( trailingOnly = TRUE )
runSeurat( args[1], args[2], args[3] )
