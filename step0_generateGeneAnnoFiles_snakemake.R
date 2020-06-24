##################################################
## Project: macaque LMD bulk RNA-seq
## Script purpose: generate gene annotation files for snakemake workflow
## Date: 2020-06-23
## Author: Yuting Liu
##################################################

## Section: set envs
##################################################
setwd("~/lustrelyt/SP/bulkRNA")
library(biomaRt)
library(refGenome)
library(plyr)

## Section: get gene, go term
##################################################
mart=useMart("ensembl")
mart <- useDataset("mmulatta_gene_ensembl", useMart("ensembl"))

ens <- ensemblGenome()
read.gtf(ens, filename = "../../monkey-brain/ref_genome/ensemble/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf")
genes <- ens@ev$genes
genes <- genes[, c("gene_name","gene_id")]
gene.ls <- unique(genes$gene_id)

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "description", 
                              "go_id", 'name_1006'),
                values=gene.ls,mart= mart)

## Section: merge go term info per genes
##################################################
df <- ddply(G_list, .(ensembl_gene_id), summarize,
             "EntreZID"=paste(unique(entrezgene_id),collapse=";"), 
             "Gene Description"=paste(unique(description),collapse=";") ,
             "GO ID"= paste(unique(go_id),collapse=";"),
            "GO Term" = paste(unique(name_1006), collapse = ";")
            )
df <- data.frame(cbind('id' = genes$gene_name[match(df$ensembl_gene_id, genes$gene_id)]),df)

# remove all the "," in the dataframe to avoid read.csv confusion
df2 <- data.frame(lapply(df, function(x) {
                  gsub(", ", " ", x)              
        }))

write.csv(df2, file = 'data/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf.annot.csv', col.names = T, row.names = F,
          quote = F)
