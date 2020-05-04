#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab
#----------------------------------------------------------------
# @Author: Mahesh Vangala (ggplot code is taken from Henry Long's original script)
# @Email: vangalamaheshh@gmail.com
# @Date: May, 9, 2016
#----------------------------------------------------------------

suppressMessages(library("reshape2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gplots"))
suppressMessages(library("RColorBrewer"))

args <- commandArgs( trailingOnly = TRUE )

LOG2TRANSFORM <- TRUE
CUTOFF <- 5

star_fusion_summary <- read.csv(args[1], header=T)
sfs <- star_fusion_summary[star_fusion_summary$TotalReads > CUTOFF,]
sfs <- sfs[,c(1,2,5)]

if (LOG2TRANSFORM) {
  sfs$TotalReads <- log2(sfs$TotalReads + 2.001)
  ds_mid = min(log2(CUTOFF + 2), max(sfs$TotalReads))
  ds_max = max(sfs$TotalReads)
  legend.title = "log2(Total\nEvidence)"
} else {
  ds_mid = min(CUTOFF, max(sfs$TotalReads))
  ds_max = sfs$TotalReads
  legend.title = "Total\nEvidence"
}

png(args[2], width = 8, height = 8, unit="in",res=300)
ggplot(suppressMessages(melt(sfs))) +
  geom_tile(aes(Sample,FusionName,fill=value), color='black') +
  scale_fill_gradient2(low = "white", high = "red",
                       midpoint = ds_mid , limit =  c(0,ds_max),
                       name = legend.title) +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0, size = 8),
        axis.text.y=element_text(size = 8)) +
  labs(title="Fusion Summary", x="",y="") +
  coord_fixed() +
  theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=1)) +
  theme(legend.position="top", legend.justification = 'right')

junk <- dev.off()
