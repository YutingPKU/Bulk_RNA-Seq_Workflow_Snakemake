suppressMessages(library("limma"))
suppressWarnings(suppressMessages(library("DESeq2")))
suppressMessages(library("edgeR"))

limma_and_deseq_f <- function(arg_counts, arg_s1, arg_s2, limma, deseq, limma_annot, deseq_annot, deseqSum_out,gene_annotation, prefix) {
    #READ in gene_annotation table--even though gene descriptions are quoted
    #in the annotations, we have to quote it again!
    if( grepl(".bz2$",gene_annotation) ) {
        gene_annot <- read.table(bzfile(gene_annotation), header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    else {
        gene_annot <- read.table(gene_annotation, header=TRUE, sep=",", quote="", fill=TRUE, row.names=NULL)
    }
    #DROP--move to the annotation files
    #Quote the gene_descriptions--VERSION 1
    #gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'],
    #                                          dQuote, simplify="vector")
    #Quote the gene_descriptions--NEED " instead of '
    gene_annot[,'Gene.Description'] <- sapply(gene_annot[,'Gene.Description'],
                                              function(x){paste0("\"",x,"\"")},
                                              simplify="vector")
    
    ## Read in lists for comparison, read in count matrix, and do "rounding" failsafe
    treatlist = strsplit(arg_s2,',')[[1]]
    ctrllist = strsplit(arg_s1,',')[[1]]
    countmat <- read.table(arg_counts, header=TRUE, sep=",", row.names=1, check.names=FALSE)
    countmat = round(countmat)
    
    ctrllist = countmat[ ,colnames(countmat) %in% ctrllist, drop = F]
    treatlist = countmat[ ,colnames(countmat) %in% treatlist, drop = F]
    
    ntreat = ncol(treatlist)
    nctrl = ncol(ctrllist)

    data = cbind(treatlist,ctrllist)
    condition = c(rep('treat',ntreat),rep('control',nctrl))
    colData <- as.data.frame(cbind(colnames(data),condition))
    
    preparationD <- function (countData, colData){
       
	dds <- DESeqDataSetFromMatrix(countData=data, colData=colData, design= ~ condition)
        dds <- dds[rowSums(counts(dds)) > 0, ]
        dds <- DESeq(dds)
        res <- results(dds)
	 
        ## Get summary stats
        summary <- c(sum(res$padj<0.05 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>1.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>1.0, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>0.0, na.rm=TRUE),
        sum(res$padj<0.05 & res$log2FoldChange < 0, na.rm=TRUE),
        sum(res$padj<0.05 & res$log2FoldChange < -0.5, na.rm=TRUE),
        sum(res$padj<0.05 & res$log2FoldChange < -1.0, na.rm=TRUE),
        sum(res$padj<0.01 & res$log2FoldChange <  0.0, na.rm=TRUE),
        sum(res$padj<0.01 & res$log2FoldChange < -0.5, na.rm=TRUE),
        sum(res$padj<0.01 & res$log2FoldChange < -1.0, na.rm=TRUE))
        sumTable <- matrix(summary, nrow=6, ncol=2)
        rownames(sumTable)<-c('log2FC > 0.0','log2FC > 0.5', 'log2FC > 1.0','log2FC < 0.0','log2FC < -0.5', 'log2FC < -1.0')
        colnames(sumTable)<-c('padj < 0.05','padj < 0.01')
        ## Write/output summary stats
        write.table(sumTable, deseqSum_out, quote=FALSE, sep=",")
        return (res)
    }
    deseq_result= preparationD(data,colData)

    ## DEG using limma; not to run
    preparationL <- function(counts,ntreat,nctrl){
        mType <- c(rep('treat',ntreat),rep('control',nctrl))
        nf <- calcNormFactors(counts)
        design <- model.matrix(~mType)
        y <- voom (counts, design, lib.size = colSums(counts)*nf)
        fit <- lmFit(y,design)
        fit <- eBayes(fit)
        res = topTable(fit,coef=2,n=length(counts[,1]),sort="p")
        return(res)
    }
    if (ntreat>1){
        limma_result = preparationL(data,ntreat,nctrl)
        limma_result <- cbind(id=rownames(limma_result), limma_result)
        ## ANNOTATE w/ local .csv biomart annotation file
        limma_annotations <- merge(limma_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
        ## MOVING to just csv files
        write.table(limma_result,limma_out,sep=',',col.names=T,row.names=F,quote=F)
        ## WRITE annotation table--limma.annot.csv
        write.table(limma_annotations,limma_annot,sep=',',col.names=T,row.names=F,quote=F)
    }
    ## Setting the first column name to 'id'
    deseq_result <- cbind(id=rownames(deseq_result), as.matrix(deseq_result))
    ## Sort by padj. and remove padj = NA
    deseq_result <-deseq_result[order(as.numeric(deseq_result[,"padj"])),]
    ## Write output deseq
    write.table(deseq_result,deseq,sep=',',col.names=T,row.names=F,quote=F)

    ## ANNOTATE w/ local .csv biomart annotation file
    deseq_annotations <- merge(deseq_result, gene_annot, by= "id", all.x = TRUE, sort=FALSE)
    
    ## defined up-regulated and down-regulated genes
    up.deseq1 <- deseq_annotations[which(deseq_annotations$log2FoldChange > 1 & deseq_annotations$padj < 0.01), c(1:9)]
    down.deseq1 <- deseq_annotations[which(deseq_annotations$log2FoldChange < -1 & deseq_annotations$padj < 0.01), c(1:9)]
    up.deseq2 <- deseq_annotations[which(deseq_annotations$log2FoldChange > 1 & deseq_annotations$padj < 0.05), c(1:9)]
    down.deseq2 <- deseq_annotations[which(deseq_annotations$log2FoldChange < -1 & deseq_annotations$padj < 0.05), c(1:9)]
    
    deseq_up1 <- paste0(prefix, ".deseq.FC2FDR001.upDEG.csv")
    deseq_up2 <- paste0(prefix, ".deseq.FC2FDR005.upDEG.csv")
    deseq_dn1 <- paste0(prefix, ".deseq.FC2FDR001.downDEG.csv")
    deseq_dn2 <- paste0(prefix, ".deseq.FC2FDR005.downDEG.csv")
    write.table(deseq_annotations,deseq_annot,sep=',',col.names=T,row.names=F,quote=F)
    write.table(up.deseq1, deseq_up1, sep = ',', col.names = T, row.names = F, quote = F)
    write.table(down.deseq1, deseq_dn1, sep = ',', col.names = T, row.names = F, quote = F)
    write.table(up.deseq2, deseq_up2, sep = ',', col.names = T, row.names = F, quote = F)
    write.table(down.deseq2, deseq_dn2, sep = ',', col.names = T, row.names = F, quote = F)
}

args <- commandArgs( trailingOnly = TRUE )
arg_counts = args[1]
arg_s1 = args[2]
arg_s2 = args[3]
limma_out=args[4]
deseq_out=args[5]
limma_annot = args[6]
deseq_annot = args[7]
deseqSum_out=args[8]
gene_annotation=args[9]
prefix=args[10]

limma_and_deseq_f(arg_counts, arg_s1, arg_s2, limma_out, deseq_out, limma_annot, deseq_annot, deseqSum_out,gene_annotation, prefix)

        
