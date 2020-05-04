# Bulk_RNA-Seq_Workflow_Snakemake
STAR, cufflinks, cluster, DESeq2, GO analysis.  
***
This workflow was designed to process bulk RNA-seq data. The output files including:

- Mapping reads with STAR;
- Counting reads with STAR and cufflinks;
- Samples quality control with PCA plot and samples-samples clustering heatmap;
- DEG with DESeq2 and limma;
- Functional enrichment analysis with GO and GSEA.

***
# metasheet.csv setting

- SampleName must be the same as config.yaml;
- Comparison group information are given by columns which start with "comp_", 1 means controls, 2 means treatments.

# config.yaml setting

- changes the samples with pathways for your input data.
Tips: generate a softlink for raw fastq files as data in the working directory.


