suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("enrichplot"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))


## The traceback is actually necessary to not break pipe at the stop step, so leave on
options(error = function() traceback(2))


goterm_analysis_f <- function(deseq_file,adjpvalcutoff,numgoterms,reference, prefix) {
    
    adjpvalcutoff = as.numeric(adjpvalcutoff)
    numgoterms = as.numeric(numgoterms)
    
    ## Read in detable
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]
 
    

    ## Append ENTREZ IDs from loaded in 
    ## Append ENTREZ IDs from loaded in db
    if (length(grep("hg",reference) == 1)) {IDdb = org.Hs.eg.db
    org = "hsa"}
    if (length(grep("mm",reference) == 1)) {IDdb = org.Mm.eg.db 
    org = "mmu"}
    detable$entrez <- mapIds(IDdb,
                             keys=rownames(detable),
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
    
    ## Select genes that pass the adjPval cutoff and select those entrez IDs as pop, set rest as universe.
    uptopgenes <- subset(detable, detable$padj < adjpvalcutoff & detable$log2FoldChange > 1)
    downtopgenes <- subset(detable, detable$padj < adjpvalcutoff & detable$log2FoldChange < -1)
    
    ## Failsafes to quit program if there are no differentially expressed genes that pass log2fc and adjpval
    if(nrow(uptopgenes) < 10) {stop(paste("Not enough significant genes for GOterm analysis, Need at least 10 and there are only ",
                                          nrow(uptopgenes)," upregulated genes at the current adjpval of ", adjpvalcutoff, sep=""))}
    if(nrow(downtopgenes) < 10) {stop(paste("Not enough significant genes for GOterm analysis, Need at least 10 and there are only ", 
                                            nrow(downtopgenes)," downregulated genes at the current adjpval of ", adjpvalcutoff, sep=""))}
    
    upselectedIDs = na.omit(uptopgenes$entrez)
    downselectedIDs = na.omit(downtopgenes$entrez)
    universeIDs = na.omit(detable$entrez)
    
    ## Perform  GO enrichment analysis
    up.bp <- enrichGO(gene = upselectedIDs, OrgDb = IDdb, ont = "BP")
    up.mf <- enrichGO(gene = upselectedIDs, OrgDb = IDdb, ont = "MF")
    up.cc <- enrichGO(gene = upselectedIDs, OrgDb = IDdb, ont = "CC")
    up.kegg <- enrichKEGG(gene = upselectedIDs, organism = org)
    
    down.bp <- enrichGO(gene = downselectedIDs, OrgDb = IDdb, ont = "BP")
    down.mf <- enrichGO(gene = downselectedIDs, OrgDb = IDdb, ont = "MF")
    down.cc <- enrichGO(gene = downselectedIDs, OrgDb = IDdb, ont = "CC")
    down.kegg <- enrichKEGG(gene = downselectedIDs, organism = org)
    
    ## visualization 
    go.ls <- c(up.bp, up.mf, up.cc, up.kegg, down.bp,down.mf, down.cc, down.kegg)
    name.ls <- apply(expand.grid("GO",c("BP","MF","CC","KEGG"),c("DEGup", "DEGdown") ), 1, paste, collapse="_")
    

    pdf(paste0(prefix,"_GOterm_dotplot_DEGFDR",adjpvalcutoff, "_TermNumber",numgoterms,".pdf"))
    mapply(function(go, name){
        g <- dotplot(go, showCategory = numgoterms) + ggtitle(name)
        capture.output(print(g))
     }, go.ls, name.ls)    
   
    dev.off()
    
    
    pdf(paste0(prefix,"_GOterm_emapplot_DEGFDR",adjpvalcutoff, "_TermNumber",numgoterms,".pdf"))
    mapply(function(go, name){
        g <- emapplot(go, title = name)
        capture.output(print(g))
    }, go.ls, name.ls)
    dev.off()
    
    up.res <- rbind(up.bp@result, up.mf@result, up.cc@result, up.kegg@result )
    down.res <- rbind(down.bp@result, down.mf@result, down.cc@result, down.kegg@result )
    
    write.table(up.res, file = paste0(prefix, "_GOterm_DEGup_results.csv"), col.names = T, row.names = F, quote = F, sep = ',')
    write.table(down.res, file = paste0(prefix, "_GOterm_DEGdown_results.csv"), col.names = T, row.names = F, quote = F, sep = ',')
    
}
gsea_analysis_f <- function(deseq_file,numgoterms,reference, prefix) {
    
    numgoterms = as.numeric(numgoterms)
    
    ## Read in detable
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]
    
    ## Append ENTREZ IDs from loaded in 
    if (length(grep("hg",reference) == 1)) {IDdb = org.Hs.eg.db}
    if (length(grep("mm",reference) == 1)) {IDdb = org.Mm.eg.db}

    
    ## prepare the genelist
    gsea_gene_list <- detable[,c("log2FoldChange")]
    names(gsea_gene_list) <- rownames(detable)
    gsea_gene_list <- sort(gsea_gene_list, decreasing=T)
    
    
    ## GSEA 
    gsea.ls <- lapply(c('BP','MF','CC'), function(vec){
        return(gseGO(geneList = gsea_gene_list, OrgDb = IDdb, ont = vec, pvalueCutoff = 0.25, keyType = 'SYMBOL'))
    })
    
    ## visualization
    name.ls <- paste0('GSEA_DEG_', c('BP','MF','CC'))
    pdf(paste0(prefix,"_GSEA_dotplot_TermNumber",numgoterms,".pdf"))
    mapply(function(go, name){
        g <- dotplot(go, showCategory = numgoterms) + ggtitle(name)
        capture.output(print(g))
    }, gsea.ls, name.ls) 
    dev.off()
    
    pdf(paste0(prefix,"_GSEA_BP_enrichplot.pdf"))
    gsea.bp <- gsea.ls[[1]]
    for (i in 1:10) {
       g <- gseaplot2(gsea.bp, geneSetID = gsea.bp@result$ID[i], title = gsea.bp@result$Description[i])
       capture.output(print(g))
    }
    dev.off()
    
    res <- rbind(gsea.ls[[1]]@result, gsea.ls[[2]]@result, gsea.ls[[3]]@result)
    
    write.table(res, file = paste0(prefix, "_GSEA_results.csv"), col.names = T, row.names = F, quote = F, sep = ',')
    
}

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
adjpvalcutoff = args[2]
numgoterms = args[3]
reference = args[4]
prefix = args[5]


goterm_analysis_f(deseq_file, adjpvalcutoff,numgoterms,reference,prefix)
gsea_analysis_f(deseq_file,numgoterms,reference,prefix)
