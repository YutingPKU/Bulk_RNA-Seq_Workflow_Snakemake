#!/usr/bin/env Rscript
# vim : syntax=r tabstop=4 expandtab 

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: May, 23, 2016
#---------------------------------

suppressMessages(library(argparse))

preprocess <- function(tpm_file, sample_names, filter_miRNA,
                       min_genes, min_samples, tpm_cutoff) {
    
	tpmTable <- read.csv(tpm_file, header=T, check.names=F,
                        row.names=1, stringsAsFactors=FALSE, dec='.')

	for (n in names(tpmTable)) {
    	tpmTable[n] <- apply(tpmTable[n], 1, as.numeric)
  	}

  	tpmTable <- na.omit(tpmTable)
  	df = tpmTable[,sample_names]
  	sub_df <- df[apply(df, 1, function(x) length(x[x>=tpm_cutoff])>min_samples),]
  	sub_df <- log2(sub_df + 1)

  	if (filter_miRNA == TRUE) {
    	sub_df <- sub_df[ !grepl("MIR|SNO",rownames(sub_df)), ]
  	}
  	min_genes = min(min_genes, nrow(sub_df))
  	## Calculate CVs for all genes (rows)
  	mean_tpm <- apply(sub_df,1,mean)
  	var_tpm <- apply(sub_df,1,var)
  	cv_tpm <- abs(var_tpm/mean_tpm)
  	## Select out the most highly variable genes into the dataframe 'Exp_data'
  	exp_data <- sub_df[order(cv_tpm,decreasing=T)[1:min_genes],]
	exp_data = exp_data[apply(exp_data, 1, var, na.rm=TRUE) != 0, ]
  	return (exp_data)
}

parse_args <- function() {
	parser <- argparse::ArgumentParser(
		description="Takes raw cuff matrix and outputs filtered one based on parameters given")
	parser$add_argument('-r', '--tpm_file', help="Raw cuff matrix file", type="character", nargs=1)
	parser$add_argument('--min_samples', help="Minimum number of samples expressed at given TPM threshold", 
				type="integer", nargs=1)
	parser$add_argument('--TPM_cutoff', help="TPM threshold", 
				type="double", nargs=1)
	parser$add_argument('--filter_miRNA', help="If TRUE, filters genes with MIR and SNO names",
				type="logical", nargs=1)
	parser$add_argument('--numgenes', help="Number of genes to be considered for plotting",
                        	type="integer", nargs=1)
	parser$add_argument('--out_file', help="Output filename to generate filtered TPM counts",
                        	type="character", nargs=1)
	parser$add_argument('--sample_names', help="Sample names from metasheet to be considered",
                        	type="character", nargs='+')
	args <- parser$parse_args()
	return (args)
}

args <- parse_args()
filtered_cuff <- preprocess(args$tpm_file, args$sample_names, args$filter_miRNA,
				args$numgenes, args$min_samples, args$TPM_cutoff)
write.csv(filtered_cuff, file=args$out_file, quote=FALSE)
