suppressMessages(library(sva))

remove_batch_effect_f <- function(countmatfile,metafile, batch_column,datatype, csvoutput,pdfoutput) {

    ## Preprocess rpkm table
    rpkmFile = countmatfile
    rpkmTable <- read.table(rpkmFile, header=T, check.names=F, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
    for (n in names(rpkmTable)) {rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)}
    rpkmTable = na.omit(rpkmTable)
    countmat = rpkmTable
    
    ## Read in metadata file
    meta <- read.table(metafile, header = TRUE, row.names = 1, sep = ",", quote = "", check.names=F)
    meta = meta[ , !grepl('comp_*', names(meta))]

    samples <- intersect(colnames(rpkmTable), rownames(meta))
    meta <- as.data.frame(meta[samples,])
    
    ## Calculate CVs for all genes (rows) - ComBat needs a minimal variance, so selecting ones that have a min var

    if (datatype == "cufflinks") {CVFILTER = 0.01}
    if (datatype == "star") {CVFILTER = 3}

    mean_rpkm_nolym <- apply(countmat,1,mean)
    var_rpkm_nolym <- apply(countmat,1,var)
    cv_rpkm_nolym <- abs(var_rpkm_nolym/mean_rpkm_nolym)
    filt_genes = subset(cv_rpkm_nolym, cv_rpkm_nolym > CVFILTER)

    ## Select those genes that pass variance filtering
    countmat = countmat[rownames(countmat) %in% names(filt_genes),]
        
    ## get known batches
    batches <- factor(meta[, batch_column])

    if(nlevels(batches) == 1) {
        # just copy the given data
        pdf(file = pdfoutput)
        dev.off()
        sink()
        countmat_nobatch <- data.frame(countmat, check.names = FALSE)
        countmat_nobatch <- cbind(Gene_ID = rownames(countmat_nobatch), countmat_nobatch)
        write.table(countmat_nobatch, file = csvoutput, sep=",", quote = FALSE, row.names = FALSE)
        print("Only one batch present, please correct metasheet")
    } else {
        # create intercept model (needed by combat, here indicating that no other covariates are involved)
        model <- model.matrix(~1, data = meta)

        # convert to matrix
        countmat <- as.matrix(countmat)

        # perform batch correction. Small values can become < 0. It is safe to max them
        # away: https://support.bioconductor.org/p/53068/
        pdf(file = pdfoutput)
        countmat_nobatch <- pmax(ComBat(countmat, batches, mod = model, prior.plots=TRUE), 0.0)
        dev.off()

        # format data.frame
        if (datatype == "star") {countmat_nobatch = round(countmat_nobatch)}
        countmat_nobatch <- data.frame(countmat_nobatch, check.names = FALSE)
        countmat_nobatch <- cbind(Gene_ID = rownames(countmat_nobatch), countmat_nobatch)

        # write to stdout
        write.table(countmat_nobatch, file = csvoutput, sep=",", quote = FALSE, row.names = FALSE)
    }
    
}


args <- commandArgs( trailingOnly = TRUE )
countmatfile = args[1]
annotFile = args[2]
batch_column = args[3]
datatype = args[4]
csvoutput = args[5]
pdfoutput = args[6]

remove_batch_effect_f(countmatfile,annotFile, batch_column,datatype, csvoutput,pdfoutput)


#remove_batch_effect_f(
#    snakemake@input[["countmat"]],
#    snakemake@input[["metadata"]],
#    snakemake@params[["batch_column"]],
#    snakemake@output[["pdfoutput"]],
#    snakemake@output[["csvoutput"]]
#)
