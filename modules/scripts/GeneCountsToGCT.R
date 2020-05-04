
##### FROM https://github.com/Bioconductor-mirror/ArrayTools/blob/master/R/output.gct.R
output.gct <- function(rpkmFile, GCTFile){
    normal = read.table(rpkmFile, sep=",", header=T, row.names=1)
    
    NAME<-Description<-rownames(normal)
    file<-cbind(NAME, Description, normal)
    #GCTFile <- paste(filename, ".gct", sep="")
    cat("#1.2", "\n", sep="\t", file=GCTFile)
    cat(nrow(normal), ncol(normal), "\n", sep="\t", file=GCTFile, append=T)
    suppressWarnings(write.table(file, file=GCTFile, row.names=FALSE, quote = FALSE, sep="\t", append=T))
}

args <- commandArgs( trailingOnly = TRUE )
rpkmFile <- args[1]
GCTFile <- args[2]

output.gct( rpkmFile, GCTFile )
