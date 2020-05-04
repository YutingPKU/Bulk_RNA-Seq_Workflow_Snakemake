#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab

#---------------------------
# @authors: Zach Herbert, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: June, 1, 2016
#---------------------------

merge_data <- function(file1, file2) {
  deseq <- read.csv(file1,header=T, row.names=1)
  limma <- read.csv(file2, header=T, row.names=1)
  result_df <- merge(deseq, limma, by="row.names")
  result_df <- result_df[,c(1,3,8)]
  return (result_df)
}

corr_plot <- function(merged_data, png_file, title="DESeq2 v limma") {
  png(png_file)
  smoothScatter(merged_data$log2FoldChange,merged_data$logFC,xlab="DESeq Log2FC", ylab="Limma LogFC",main=title)
  abline(lm(merged_data$log2FoldChange ~ merged_data$logFC),col="red")
  abline(v=0,h=0,lty=2)
  junk <- dev.off()
}

args <- commandArgs(trailingOnly=TRUE)
result_df <- merge_data(args[1], args[2])
write.csv(result_df, file=args[3])
corr_plot(result_df, args[4])
