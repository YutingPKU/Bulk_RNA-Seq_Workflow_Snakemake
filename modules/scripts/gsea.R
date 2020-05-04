#Script to run enricher fn on most significantly diffexp genes

suppressMessages(library("clusterProfiler"))
suppressMessages(library(ggplot2))

gsea <- function(deseqTable, gsea_db, comp_title, out_path) {
   #Constant used to determine how many of the top hits to use
   N <- 10

   #get set of significant diffexp genes--threshold 0.05
   genes <- deseqTable[deseqTable$padj < 0.05,]

   #sort by padj --DROP b/c this is already sorted by padj
   #genes <- genes[order(genes$padj),]

   #get list of gene names
   gene_list <- rownames(genes)
   
   #WRITE this as out_gene_list
   #write.list(gene_list, "gene_list.txt", quote=F, col.names = NA, sep="\t")
   write(gene_list, file=paste0(out_path, ".gene_list.txt"), ncolumns=1)

   #CHECK if the egmt output already exists
   out_egmt = paste0(out_path, ".egmt.Rda")
   if (file.exists(out_egmt)) {
      #print("Exists")
      #LOAD
      load(out_egmt)
   } else {
      #print("New")
      #START from scratch
      #running enrichment on gene symbol
      #READ in gmt
      gmt <- read.gmt(gsea_db)
      egmt <- enricher(gene_list, TERM2GENE=gmt)
      save(egmt, file=out_egmt)
   }
   #WRITE out summary
   write.csv(summary(egmt), file=paste0(out_path, ".gene_set.enrichment.txt"))

   #generate dotplot
   png(paste0(out_path, ".gene_set.enrichment.dotplot.png"), width = 8, height = 8, unit="in",res=300)
   #NOTE: if I just did dotplot(...), the png doesn't get saved.  
   #I need the print(p)
   
   topResults = egmt@result[1:N,]
   #NOTE: the $GeneRatio is an character array, e.g. "120/1317" we want to
   #convert those to float.  This might be non-canonical, but we're simply
   #evaluating them and saving them to a new col, GeneRatios 
   #NOTE: GeneRatio = original char expr, GeneRatios = evaluated!
   topResults$GeneRatios <- apply(as.array(topResults$GeneRatio), 1, function(expr) eval(parse(text=expr)))

   p <- ggplot(topResults, aes(x = reorder(topResults$Description, -log(topResults$qvalue)), y = -log(topResults$qvalue))) +
           geom_point(aes(size = Count)) + theme_bw(base_size=10) + coord_flip() +
           labs(y="-Log(qvalue)", x="Gene Sets") +
           ggtitle(comp_title)
   print(p)
   junk <- dev.off()
   
   #---------------------------------------------------------------------------
   #GSEA plot- generate gsea list enrichment plots for top 10 gene lists
   #---------------------------------------------------------------------------
   #CHECK if the egmt output already exists
   out_gsea = paste0(out_path, ".gsea.Rda")
   if (file.exists(out_gsea)) {
      #LOAD
      load(out_gsea)
   } else {
      #running GSEA on gene_list
      #CHECK if gmt is already read-in and in memory-
      if (!exists("gmt")) {
         #READ in gmt
      	 gmt <- read.gmt(gsea_db)
      }
      #MAKE the ordered gene list
      gsea_gene_list <- genes[,c("log2FoldChange")]
      names(gsea_gene_list) <- rownames(genes)
      #SORT gsea_gene_list in DESCENDING padj --i.e. least sig. first
      gsea_gene_list <- sort(gsea_gene_list, decreasing=T)

      #make TERM2GENE mapping
      egmtTERM2GENE=gmt[,c("ont","gene")]
      
      #GSEA call
      gsea_result <- GSEA(gsea_gene_list, TERM2GENE=egmtTERM2GENE)
      save(gsea_result, file=out_gsea)
   }
   
   #for the top 10 hits, generate a (gene set) enrichment plot
   top_hits<- rownames(summary(egmt)[1:10,])
   for (gene_set in top_hits) {
      #print(gene_set)
      png(paste0(out_path, ".gsea.", gene_set,"%02d.png"), width = 8, height = 8, unit="in",res=300)
      p2<-gseaplot(gsea_result, gene_set)
      print(p2)
      dev.off()
   }
}

## Read in arguments
args <- commandArgs( trailingOnly = TRUE )
deseqFile=args[1]
gsea_db=args[2]
title=args[3]
out_path=args[4]

deseqTable <- read.csv(deseqFile, header=T, check.names=F, row.names=1, stringsAsFactors=FALSE, dec='.')
gsea(deseqTable, gsea_db, title, out_path) 
