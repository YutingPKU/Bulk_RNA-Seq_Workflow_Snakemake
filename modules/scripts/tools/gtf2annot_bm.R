#read in the file args
args <- commandArgs( trailingOnly = TRUE )
#arg_gtf = args[1]
arg_dataset = args[1]
#arg_dataset = "hsapiens_gene_ensembl"
arg_out = args[2]
if (length(args) != 2) {
    print("USAGE: Rscript gtf2annot.R [biomart dataset, e.g. hsapiens_gene_ensembl or mmusculus_gene_ensembl] [output filename]");
    quit();
}

library(biomaRt)

#ADD addition fields by biomart
mymart = useMart(biomart="ensembl", dataset=arg_dataset);
gene_annots <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "external_gene_name", "description", "go_id", "name_1006"), mart=mymart)
colnames(gene_annots) <- c("id", "EnsemblID", "EntrezID", "Gene Description", "GO ID", "GO Term")
print(length(gene_annots))
print(head(gene_annots))
print(head(gene_annots$id))
print(length(gene_annots$id))
print(length(unique(gene_annots$id)))
gene_annots[,'Gene Description'] <- sapply(gene_annots[,'Gene Description'],
                                           function(x){paste0("\"",x,"\"")},
                                           simplify="vector")   
#CHECK for NA??
write.table(gene_annots,arg_out,sep=',',col.names=T,row.names=F,quote=F)
