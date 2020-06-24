# MVIPER: Bulk_RNA-Seq_Workflow_Snakemake
STAR, cufflinks, cluster, DESeq2, GO analysis.  
***
This workflow was designed to process and visualize bulk RNA-seq data. The output files including:

- Mapping reads with STAR;
- Counting reads with STAR and cufflinks;
- Samples quality control with PCA plot and samples-samples clustering heatmap;
- DEG with DESeq2 and limma;
- Functional enrichment analysis with GO and GSEA.

***
# Table of Contents
1. [MVIPER](#mviper)  
2. [Working directory structure](#working-directory-structure)  
3. [How to run the MVIPER](#how-to-run-the-mviper)  
4. [Running VIPER](#running-viper)
5. [Outputs of MVIPER](#outputs-of-mviper)


## MVIPER
MVIPER is a bulk RNA-seq analysis pipeline built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home). MVIPER is modified [VIPER](https://bitbucket.org/cfce/viper/src/master/). Modifications are as the follows:
- add scripts to generate config.yaml and metashee.csv automatically
- add script to generate gene annotation file which is used in the ref.yaml for new species
- modified *file_format.snakefile* to transfer bam file to BigWig file format using [deeptools](https://deeptools.readthedocs.io/en/develop/)
- modified *preprocess.snakefile* to set multiple thresholds for top variant genes lists
- modified *DE.snakefile* to set multiple thresholds for DEG selection
- modified serveral R scripts and python scripts to get more beautified plots

Cornwell M, Vangala M, Taing L, Herbert Z, Köster J, Li B, Sun H, Li T, Zhang J, Qiu X, Pun M, Jeselsohn R, Brown M, Liu XS, Long HW. VIPER: Visualization Pipeline for RNA-seq, a Snakemake workflow for efficient and complete RNA-seq analysis. BMC Bioinformatics. 2018 Apr 12; 19(1):135. PMID: [29649993](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/pubmed/?term=29649993).

## Working directory structure
Download the MVIPER by the following command:
`git clone https://github.com/YutingPKU/Bulk_RNA-Seq_Workflow_Snakemake.git`

Put the MVIPER and your input data in a directory (like PROJECT).
> PROJECT/   - *the root directory*
> modules/   - *the scripts and snakefile for MVIPER*
> static/ - *the reference metadata for MVIPER*
> data/   - *the input data directory*
> config.yaml   - *pipeline configure file*
> ref.yaml - *reference metadata configure file*
> metasheet.csv - *metadata for input data samples*
__Note__: the *data* directory is the pathway of your input data, you can generate a soft link to your fastq files. Example is listed in the following
`ln -s fastq_directory data`

## How to run the MVIPER
1. __setting the ref.yaml__
-  star_index: the genome index directory for [STAR](https://github.com/alexdobin/STAR) alignment
-  gtf_file: gene annotation file by gtf file format
- gene_annotation: csv file containing gene symbol, gene description, ENSEMBL id, EntreZ id, GO id and GO term. For human, mouse and macaque, you can use the bz2 file from the *static/* directory. For other species, you can generate the gene annotation file by the following command:
`Rscript step0_generateGeneAnnoFiles_snakemake.Rstep0_generateGeneAnnoFiles_snakemake.R `
__Note__: change the biomart dataset according to your species
2. __setting the config.yaml__
- set the configure to control the pipeline running and scripts parameter
- adding the  sample names and pathways(fastq file directory) of the input data to the *samples* key. Examples are listed in the following
```
 samples:
    10068A-CPi-RNA-lib:
      - data/10068A-CPi-RNA-lib/10068A-CPi-RNA-lib_R1.fq.gz
      - data/10068A-CPi-RNA-lib/10068A-CPi-RNA-lib_R2.fq.gz
```
You can add the *samples* automatically by the following command
    `bash addDataInfo_configymal.sh`

3. __setting the metasheet.csv__
- SampleName: sample names of the input data, must be the same as config.yaml
- You can add any annotation information for the samples by adding in the columns
- Comparison group information are given by columns which start with "comp_", 1 means controls, 2 means treatments. You can set multiple comparison groups. Examples are listed in the following

| SampleName | Regions | Replicates | comp_CPivsOther | comp_CPovsOther |
|--------|------|------------|-----------|------------|---------------|---------------|
|10068A-CPi-RNA-lib   | CPi | 10068A | 2     |   1          |
|10068A-CPo-RNA-lib   | CPo | 10068A | 1    |   2          |
4. __test and run the workflow__
- validate the pipeline by
`snakemake --np -s viper.snakefile`
- run the pipeline by
`snakemek -s viper.snakefile --cores 20 -j 10`
- run the pipeline in cluster by
`snakemake -j 10  -pr  -c "pkubatch -p cn-short -N 1 -c 20 --qos=lch3000cns -A lch3000_g1   -J {rule}.{wildcards} -o logs/cluster/{rule}/{rule}.{wildcards}_%j.out -e logs/cluster/{rule}/{rule}.{wildcards}_%j.err " -s viper.snakefile -k 2`
__Note__: setting the account information according to your cluster account

##Outputs of MVIPER
benchmarks files and log files are in the *bednchmarks* and *logs* directories. Analysis results and Visualization resutls are in the *analysis* directory. The Contents of *analyais* directory are listed in the following
```
bam2bw/ - BigWig files of the aligned bam files per sample
cufflinks/ - genes and isoforms fpkm matrix per sample
STAR/ - STAR alignment results per sample
summary_reports/
  cufflinks/ - genes fpkm matrix for all samples; top variant genes fpkm matrix for all samples
  STAR/ - alignment reports for all samples
  plots/ - clustering heatmap based on samples-samples correlation matrix; PCA plots based on top variant genes fpkm matrix
  diffexp/ - DEG lists detected by DESeq2 and limma; Vocalno plots for DEG
  GO/ - GO, KEGG Pathway and GSEA analyses and Visualization for DEG
```

## Contact
If you have any questions about MVIPER, please feel free to contact lyt17@pku.edu.cn
