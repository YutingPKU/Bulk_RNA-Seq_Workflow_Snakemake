#!/usr/bin/env Rscript

# vim: syntax=r tabstop=4 expandtab

#---------------------------------
# @authors: Zach Herbert, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: Aug, 05, 2016
#---------------------------------

library(reshape2)
library(ggplot2)

options(error = function() traceback(2))
args <- commandArgs( trailingOnly = TRUE )

cuff.csv <- args[1]
plot.png <- args[2]

cdat <- read.csv(cuff.csv, header=T, row.names=1)
fpkm_above_0.1 <- colSums(apply(cdat, 2,function(x) ifelse(x>0.1, 1, 0)))
fpkm_above_1 <- colSums(apply(cdat, 2,function(x) ifelse(x>1, 1, 0)))

fpkm <- cbind(fpkm_above_0.1,fpkm_above_1)
mfpkm <- melt(fpkm)

yl <- 1.2*(max(mfpkm$value))


png( plot.png,width = 8, height = 12, unit='in',res=300)
ggplot(mfpkm,aes(Var1,value,fill=Var2))+
  ggtitle("\nGenes Detected\n")+
  xlab("Samples\n")+
  ylab("\nNumber of Genes")+
  geom_bar(stat="identity",
           position="dodge",
           width=.6,color="white")+
  scale_fill_manual(values=c("slategray3","plum3"))+
  geom_text(aes(label=value),
            vjust=-0.2, 
            position=position_dodge(width=0.5),
            hjust=-0.25, size=4)+
  coord_cartesian(ylim=c(0,yl))+
  #facet_grid(Var2~.)+
  theme_bw()+
  theme(text = element_text(size=12), 
        axis.text.x = element_text(face = "bold",angle=0,size=16),
        axis.text.y = element_text(face = "bold",size = 12),
        axis.title = element_text(color="black",size = 18, face = "bold"),
        plot.title = element_text(lineheight=.8, size = 20, face="bold"),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(colour="black", size = 16, face = "bold"),
        #strip.text.y = element_text(color="black",size = 16, face = "bold"),
        #strip.background = element_rect(colour="black",fill="lightgrey",size=1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
        ) +
  coord_flip()

junk <- dev.off()


