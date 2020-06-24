
# VIPER - Visualization Pipeline for RNA-seq

# Introduction to VIPER:

__VIPER__ is a comprehensive RNA-seq analysis tool built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users. VIPER combines the use of several dozen RNA-seq tools, suites, and packages to create a complete pipeline that takes RNA-seq analysis from raw sequencing data all the way through alignment, quality control, unsupervised analyses, differential expression, and downstream pathway analysis. In addition, VIPER has been outfitted with several recently published tools that allow for interrogation of immune and virus infiltrate. The results are compiled in a simple and highly visual report containing the key figures to explain the analysis, and then compiles all of the relevant files, tables, and pictures into an easy to navigate folder.

Cornwell M, Vangala M, Taing L, Herbert Z, Köster J, Li B, Sun H, Li T, Zhang J, Qiu X, Pun M, Jeselsohn R, Brown M, Liu XS, Long HW. VIPER: Visualization Pipeline for RNA-seq, a Snakemake workflow for efficient and complete RNA-seq analysis. BMC Bioinformatics. 2018 Apr 12; 19(1):135. PMID: [29649993](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/pubmed/?term=29649993).

# Table of Contents
0. [System Requirements](#requirements)
1. [Anatomy of a VIPER PROJECT](#anatomy)  
2. [Getting Started - One time Installation of components necessary for an individual user](#GettingStarted)  
3. [Setting up a Project Folder for a VIPER Run](#SettingUpForProject)  
	1. [The Config File](#config)  
	2. [The Metasheet](#metasheet)  
4. [Running VIPER](#RunningViper)
5. [Appendix A: Dana-Farber Members](#DFmembers)
6. [Appendix B: Specific VIPER Commands for Replotting](#replotting)
7. [Setting up VIPER for a group of users or server](#serverSetup)

# System requirements:
Some of the tools that VIPER uses, e.g. STAR and cufflinks are very memory intensive programs.  Therefore we recommend the following system requirements for VIPER:

### Minimal system requirements:
We recommend that you run VIPER on a server that has at least 30GB of ram.  This will allow for a single-threaded VIPER run (on human samples).

### Recommended system requirements:
We recommend that you have at least 128GB of ram and at least a 4-core CPU if you want to run VIPER in multi-threaded mode (which will speedup the workflow significantly).  Our own servers have 256GB of ram and 32 cores.

# Anatomy of a VIPER PROJECT: <a name="anatomy"></a>
All work in __VIPER__ is done in a __PROJECT__ directory, which is simply a directory to contain a single __VIPER__ analysis run.  __PROJECT__ directories can be named anything (and they usually start with a simple mkdir command, e.g. mkdir viper_for_thesis),  but what is CRITICAL about a __PROJECT__ directory is that you fill them with the following core components:
(We first lay out the directory structure and explain each element below)
> PROJECT/  
> viper/  
> data/  - *optional*   
> config.yaml  
> metasheet.csv

The 'viper' directory contains all of the viper code.  We'll explain how to download that directory below.  The 'data' directory is an optional directory that contains all of your raw data.  It is optional because those paths __may__ be fully established in the config.yaml, __however__ it is best practice to gather your raw data within 'data' using [symbolic links](https://www.cyberciti.biz/faq/creating-soft-link-or-symbolic-link/).

The *config.yaml* and *metasheet.csv* are configurations for your VIPER run (also explained below).

After a successful __VIPER__ run, another 'analysis' folder is generated which contains all of the resulting output files.

# Getting Started - One time installation of components necessary for an individual user: <a name="GettingStarted"></a>
__If you are looking to install for a system of users, we recommend you look at appendix C below. Note that this can also be a very useful step for individual users as well!__

Although included in this README are step-by-step instructions, it is assumed that the user has a basic understanding of the [nix command line interface](https://en.wikipedia.org/wiki/Command-line_interface).

### Installing wget and git:

To get some of the required software packages, we will use the command line tools called [wget](http://www.gnu.org/software/wget/) and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).  *wget* is a popular tool for downloading things off of the internet.  *git* is a distributed version control system which we will use to checkout the VIPER code.

__These tools are already pre-installed in most systems__, but if you are unsure whether or not you have *wget* enter `wget` and if the return is `wget: command not found`, then you will have to install *wget*.  Do likewise for *git*.

### Installing Miniconda3:

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system (aka __CONDA__) to manage all of the software packages that __VIPER__ is dependent on. 

Use following commands to retrieve and then __RUN__ the Minicoda3 installation script:  
1.	`wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
2.	`bash Miniconda3-latest-Linux-x86_64.sh`  

- Whilst running the installation script, follow the commands listed on screen, and press the _enter_ key to scroll.
- __Make sure to answer yes when asked if you want to prepend Miniconda3 to PATH.__
- Close your terminal, open a new one and you should now have Conda working! Test by entering:  
	`conda update conda`
	- Press `y` to confirm the conda updates

__NOTE__: you will only have to install Minicoda3 once.  
__NOTE__: remember to close your terminal session and re-login

### Installing the VIPER conda environments:

We are now ready to use __CONDA__ to install the software packages which __VIPER__ is dependent on.

1.	`wget https://bitbucket.org/cfce/viper/get/master.tar.gz`
2.	`tar -xf master.tar.gz`
3.	`mv cfce-viper-XXXXX viper`  
__NOTE__: the XXXXX refers to the latest changeset of viper, so it will differ  
4.	`cd viper/envs`
5.	`conda env create -f environment.yml -n viper`
6.	`conda env create -f python2_environment.yml -n viper_py2`

__NOTE__: you will only have to install the VIPER conda environments once.

### DOWNLOADING the VIPER static reference files:

__VIPER__ is dependent on reference files which can be found for the supported species listed below:  [download link](https://www.dropbox.com/sh/8cqooj05i7rnyou/AAB-i4hHxQwqJDTXbzM_2JPua?dl=0)
hg19
mm9

To unzip these files: 
`tar -xzf hg19.tar.gz` OR `tar -xzf mm9.tar.gz`

__BEST PRACTICE:__ we recommend that you download the reference files that you need and then untarring then in a directory called "VIPER_static".  So for example, suppose you make "VIPER_static" in you home directory then you would have the following directory structure:
> VIPER_static/  
>		hg19/  
>		mm9/  

__NOTE__: you will only have to download the static references once.

# Setting up your PROJECT folder for a VIPER run: <a name="SettingUpForProject"></a>

We are now ready to setup __VIPER__ within our __PROJECT__ folder.  
If you haven't already, start by making your __PROJECT__ folder and changing into that directory.  

`mkdir PROJECT`  
`cd PROJECT`  

__NOTE__: remember, you can name your __PROJECT__ folder anything

### Linking the VIPER static reference files (optional best practice):

If you followed the __BEST PRACTICE__ of placing the reference files in a 
directory called "VIPER_static" then you __should__ create a symbolic link to that directory and name it 'ref_files':  
1.	`cd PROJECT`  
2.	`ln -s /path/to/VIPER_static ./ref_files`  
__NOTE__: if VIPER_static is in your home directory, then the path would be ~/VIPER_static

This is the recommended __BEST PRACTICE__ because it will ease setting up the __PARAMETERS__ and __PATHS__ in config.yaml described below.  If you don't choose to implement this __BEST PRACTICE__ then you can manually define the paths to the required files in config.yaml.

### DOWNLOADING the VIPER source code:
Within your __PROJECT__ directory, issue the following commands:  
1.	`wget https://bitbucket.org/cfce/viper/get/master.tar.gz`  
2.	`tar -xf master.tar.gz`  
3.	`mv cfce-viper-XXXXX viper`  
__NOTE__: the XXXXX refers to the latest changeset of viper, so it will differ  

__ADVANCED__: you may clone the latest version of [__VIPER__](https://bitbucket.org/cfce/viper) using git

### Setting up your 'data' directory (optional):

As mentioned above, we highly recommend that you pool your raw data into a directory called 'data' (within __PROJECT__).  

__IF all of your data is already centrally stored__ in a directory called '/some/path/to/my/data/', then symbollically link your data to your PROJECTS folder.  
`cd PROJECT`  
`ln -s /some/path/to/my/data/ ./data`  
The 'ln -s' command creates a symbolic link from '/some/path/to/my/data/' and names it 'data'. Symbolic Links have been shown to be adequate when testing, if not, then please try one of the other solutions.

__IF your files are not centrally stored__:  
- make a directory within __PROJECT__ called 'data' and copy your raw files into data

### Copying over the __META__ files:

The __META__ files (*config.yaml* and *metasheet.csv*) allow you to configure each run.  They are explained in much more detail below.  For now, we simply copy them from the viper source directory:
```
	cd PROJECT
	cp viper/config.yaml .
	cp viper/metasheet.csv .
```
__We will explain how to edit and configure these files shortly below__

##### What your PROJECT directory should look like (up to now):
> PROJECT/  
> ref_files/ - optional, but __BEST_PRACTICE__ is highly recommended*  
> viper/  
> data/  - *optional*  
> config.yaml  
> metasheet.csv  

__NOTE__: you will have to setup your PROJECT directory for each VIPER run.

## Configuring the META files: config.yaml <a name="config"></a>
The config.yaml file has three main sections. __PATHS__, __PARAMS__, __SAMPLES__:

#### PATHS:
In this section, you will need to specify the location of the following static reference files.  

__BEST PRACTICE:__ IF you followed the best practices of downloading the VIPER static reference files to a directory in your home directory named "VIPER_static" and then created a symbolic link to this directory (calling the link "ref_files") then you can __skip most__ of this section __As the default config.yaml ASSUMES that ref_file points contains the VIPER STATIC reference directories__ (i.e. anything that starts with ./ref_files can be skipped/ignored).  Simply make sure that your references are appropriate for your species/assembly.  The only parts you will need to be responsible for are downloading the ref_fasta file and generating the STAR index as these were too large to include in the static reference files. 

__NOTE:__ Even if you haven't followed the __BEST PRACTICE:__, once you've configured your paths for your *first* run you can usually just use it as a template for future runs by copying it. __SO this is something you probably will only have to do once!__

>bed_file: ./ref_files/hg19/RefGene/refseqGenesHg19.bed  
>  -  RefSeq gene file for your assembly
>
>genome_lib_dir: ./ref_files/hg19/Hg19_CTAT_resource_lib  
>  -  path to your CTAT_resource_lib (used by STAR-Fusion module)
>
>gtf_file: ./ref_files/hg19/Hg19_CTAT_resource_lib/ref_annot.gtf  
>  -  Gene annotation file
>
>
>ref_fasta: /some/path/to/humanhg19/rawgenome/hg19.fasta  
>  -  Reference genome in a .fasta or .fa format
>  -  __IMPORTANT__ Assemblies are not found within the VIPER_static reference files as they are prohibitively __LARGE__.  You can download genome assemblies here at [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html)
>  -  __IMPORTANT__ You must also index the fasta file by issuing this command  
>  -  `samtools faidx [fasta file]`, e.g. samtools faidx hg19.fasta
>
>reference: hg19  
>  -  What assembly (shortname) you are using, current options are hg19 and mm9
>
>star_index: /some/path/to/STAR/humanhg19/  
>  -  star-index for the STAR aligner  
>  -  __IMPORTANT__ STAR indices are not found within the VIPER_static reference files as they are prohibitively __LARGE__.  
>  -  You can generate one by running the following:  
>  -  `STAR  --runMode genomeGenerate --runThreadN 24 --genomeDir /where/you/want/to_store/the/output -genomeFastaFiles /dir/to/hg19/hg19.fa`  
>    - __runThreadN__ 24 means to use 24 cores (optional parameter)  
>
>star_rRNA_index: ./ref_files/hg19/humanhg38_ncrna/  
> - Follow same instructions as for normal star-index
> 
> gene_annotation: viper/static/humanhg19.annot.csv  
>  -  Path to annotation files (e.g. ENSEMBL\_ID, Gene Description, Go Terms, etc.).  This file can be generated by using [SCRIPT that Len needs to include].  Pre-made annotations for hg19 and mm9 can be found in viper/static (simply bunzip2 them).
>  -  __Simply select the appropriate .annot.csv in viper/static appropriate for your species__


#### PARAMS:

This section holds parameters specific to your project design

>metasheet: metasheet.csv  
>  - The name and path of your metasheet. This should be within the PROJECT directory as described above. If you decide to change the name of your metasheet so it is specific to your project, you will have to change it here. for example, it is good practive to name it something like mutation\_metasheet.csv
>
>stranded: true  
>  -  MISSING
>    
>library_type: 'fr-firststrand'  
>  -  MISSING
>
>RPKM\_threshold: 1.0  
>  -  Minimum RPKM Value for a gene to be significant
>
>min\_num\_samples\_expressing\_at\_threshold: 4  
>  -  Number of samples that need to have the minimum RPKM threshold for the gene to be significant
>
>filter\_mirna: true  
>  -  Filter out MicroRNA (RNA that start with "SNO" or "MIR")
>
>numgenes_plots: 1000  
>  -  Number of Genes to be shown in the plots of VIPER, including PCA, Sample-Sample, and Sample-Feature
>
>num\_kmeans\_clust: 0,4  
>  -  Parameter for the Sample-Feature Clustering Heatmap. VIPER will interpret this list of numbers as the types of SF heatmaps the user wants. 0 means hierarchical, 4 means kmean of 4, 0,4 means both a hierarchical and a 4 kmeans clustered heatmap
> 
> snp_scan_genome: false  
>  -  Boolean Flag {true | false} on whether to perform a genome-wide snp scan *IN ADDITION* to the snp scan done on chr6.
>
> cancer_type: "sarc"  
> - Tells VIPER to perform a TIMER analysis which will generate an estimate of immune cell type abundances in your cancer samples.  
> - __NOTE__: Your samples must be one of the following TCGA cancer types:  
> - Cancer types available {'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg', 'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct', 'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca', 'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'}  
> - See https://cistrome.shinyapps.io/timer/ ('Cancer Type') drop-down for type descriptions.
> - To ENABLE: uncomment and put in your sample's TCGA cancer type.  
> - To DISABLE: Comment out if not needed or set to 'False'! (default)  
>
> cdr3_analysis: false  
> - Boolean Flag {true | false} Tells VIPER to TRUST analysis which try to determine the CDR3 sequences in your __PAIRED END__ samples.  
> - To ENABLE: uncomment and set to true  
> - To DISABLE: Comment out if not needed or set to false! (default)  
>
> virus_dna_scan: false  
> - Boolean Flag {true | false} Tells VIPER to detect viral dna sequences in your human samples.  
> - To ENABLE: uncomment and set to true  
> - To DISABLE: Comment out if not needed or set to false! (default)  


##### SAMPLES:

In this section of the configuration file, you specify the __NAMES__ of each sample, and the __PATHS__ to the sample's raw data.  Raw data files can either be fastq, fastq.gz, or bam formated files.

As recommended above, if all of your raw data are located in __PROJECTS/data__, then each path will simply start like:  
`'data/first.fastq'`

__If you did not follow the recommended best practice__ then you will have to specify the full paths here.

Each sample should be given a __NAME__ (arbitrary text) and a __PATH__

__EXAMPLE__:
```
samples:
	SAMPLE1:
		- data/SAMPLE1.fastq.gz
	SAMPLE2:
		- data/SAMPLE2.fastq.gz
```  

###### __For Paired-end samples, simply add the second samples of the pait__

__EXAMPLE__:
```
samples:
	SAMPLE1:
		- data/SAMPLE1_R1.fastq.gz
		- data/SAMPLE1_R2.fastq.gz
	SAMPLE2:
		- data/SAMPLE2_R1.fastq.gz
		- data/SAMPLE2_R2.fastq.gz
```

__IMPORTANT__: __You cannot mix Paired-end and Single-end samples within the same VIPER run as this will cause an ERROR__. If necessary, run all of one type of data first, followed by the net type of data after.


## Configuring the META files: metasheet.csv <a name="metasheet"></a>

Make the *__metasheet__* file in excel, and save it as a .txt or .csv, It doesn’t matter what it is named as long as it is called in the *__config__* in the spot marked “metasheet,” see the *__config__* section if confused. The format should be something like the following:

| Sample | Cell | Condition  | Treatment | Replicates | comp_MCF7_AvB | comp_T47D_CvD |
|--------|------|------------|-----------|------------|---------------|---------------|
| A1     | MCF7 | Full_Media | NoDOX     | 1          | 1             |               |
| A2     | MCF7 | Full_Media | NoDOX     | 2          | 1             |               |
| B1     | MCF7 | Full_Media | DOX       | 1          | 2             |               |
| B2     | MCF7 | Full_Media | DOX       | 2          | 2             |               |
| C1     | T47D | Full_Media | NoDOX     | 1          |               | 1             |
| C2     | T47D | Full_Media | NoDOX     | 2          |               | 1             |
| D1     | T47D | Full_Media | DOX       | 1          |               | 2             |
| D2     | T47D | Full_Media | DOX       | 2          |               | 2             |

- The first column should always be sample names that exactly match the sample names used in config.yaml (see __SAMPLES__ just above)
- The samples that you want to perform a Differential Expression (DE) on using limma and deseq should be marked by the “comp” columns more on this below
	- This is important! The “control” should be marked with a 1, and the “treatment” should be marked with a 2.
- It is recommended that if you should have a “replicates” column to denote different samples, it is a good idea to not only have each of the sample names be unique, but also make sure that the associated metadata is unique as well to each sample.
- The rest of the  metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many “comp_XXXX” columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

- Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.
- Common Problems with *__metasheet__*
- Characters to avoid: ("-", "(", ")", " ", "/", "$") To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation, or they are invalid characters in R. NOTE: viper parses the meta file and will convert MOST of these invalid characters into '.'--dollarsigns will just be dropped.  The viper parser will also convert between dos/mac files to unix format.
	- It is very important that you know that samples A is what you mark with 1, and samples B is what you mark with a 2. You should name your output following this format as well "comp\_cond\_AvB” This will let the reader know what the output DE files refer to.  
		-  Deseq: ”baseMeanA” refers to samples A, which follows condition 1 and “baseMeanB” refers to samples B which follows condition 2. logfc is B/A
		-  Limma: Logfc refers to B/A

# Running VIPER: <a name="RunningViper"></a>

Now that we have setup our __PROJECT__ directory (downloading the 'viper' code directory, creating our 'data' directory, and configuring our config.yaml and metasheet.csv), __we are (finally!) ready to run VIPER.__

To start, we must activate the __VIPER CONDA ENVIRONMENT__:
`source activate viper`  
*if successful, you will see "(viper)" prepended to your command prompt.  
__NOTE__: you will have to do this for every VIPER run/session.  Once you log-out of your terminal session, you will also log out of your VIPER CONDA ENVIRONMENT.

Next we will perform a __DRY-RUN__ to make sure that we setup the VIPER PROJECT directory correctly.  In your __PROJECT__ folder run the following command:

`snakemake --snakefile viper/viper.snakefile -n`  

This will return a large output which basically outlines what VIPER is about to do. If no errors come back, then you will mostly see __GREEN__ and __YELLOW__ print outs.  If there are errors, then you will see some __RED__ print outs.

__If there are no errors, then use the following command to run VIPER:__

`snakemake --snakefile viper/viper.snakefile`

*If there are errors, try to see what the error is about.  Was it a mistyped path?  Etc.  If all else fails, email the VIPER team (email address needed)*

### APPENDIX A: Dana-Farber CFCE Members: <a name="DFmembers"></a>
If you are a member of Dana-Farber and have access to the CFCE server, you will already have many of the packages you need installed globally. Please see the [README within the CFCE folder of VIPER](https://bitbucket.org/cfce/viper/src/772915c62ff08ff951f813746d277fbe60f71a45/cfce/README_CFCE.md?at=master&fileviewer=file-view-default)


### APPENDIX B: Specific Replotting: <a name="replotting"></a>
After you have run __VIPER__ in its entirety, you may want to go back and tweak your outputs. Maybe adding or subtracting metadata columns, differential expression columns, or maybe just doing a subset of your data. Below is a list of snakemake commands to run __VIPER__ to rerun some specifics for further downstream analysis tweaking. Note that __VIPER__ is built to automatically rerun downstream analysis if you adjust the *config* or the *metasheet*.

To learn about how snakemake works, and some of the specifics of the following commands and others, look into the [snakemake documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation)

The following are some useful commands for rerunning and adding to the download analysis without having to rerun the whole pipeline:

`snakemake -s viper/viper.snakefile -n`  

`snakemake -s viper/viper.snakefile -j 24`

`snakemake -s viper/viper.snakefile analysis/plots/heatmapSF_plot.pdf -f `
  
`snakemake -s viper/viper.snakefile analysis/plots/heatmapSS_plot.pdf -f `

`snakemake -s viper/viper.snakefile analysis/plots/pca_plot.pdf -f `

Adding comp columns will automatically make it generate new differential expressions analysis and adjust figures accordingly.  
"Touching" the metasheet will have __VIPER__ rerun all downstream output. 

`touch metasheet.csv`  

### APPENDIX C: Deploying for a group of users: <a name="serverSetup"></a>
NOTE: this section is by no means "the solution".  It is just the particular solution that we deployed for our center.

Problem1: The install instructions above applied to single users.  But suppose you work within a lab and you want all of your lab users to use VIPER.  You want a way to centralize your VIPER deployment so that everyone in your lab is using the same (updated) version of VIPER.

Problem2: Installing miniconda for in your own local environment has a tendency to "clobber" the global enviroment on your system.  For example, your lab/server uses python2.7, but once your install the latest miniconda pkgs, you find yourself *forced* to use python3.5.  Or on the server, STAR ver X.Y is installed, but with conda, you are using STAR ver Z.A.

There are ways to hack around this by modifying your .bashrc but you don't even know what a .bashrc is.

Solution: One solution to these problems is to install VIPER (as above) for one "central" user--we created a 'viper' user account on our ubuntu system. And then we created a simple initialization script, 'viperSetup.sh', which setup the PROJECT directory for our users.

Let's walk through this setup in more detail:

### setting up a central viper user:
1. on your machine, create a new user--e.g. 'viper'
2. check out the latest VIPER code: e.g.
   /home/viper$ `git clone git@bitbucket.org:cfce/viper.git`
3. install miniconda and create the VIPER conda enviromentS, etc: 
   see "Installing VIPER and setting up your environment" above
4. optional but recommended: create template config.yaml and metasheet.csv that will help your users run VIPER.  For example, in this template config.yaml, preset all of the paths that VIPER will need or comment out/in the options that will be commonly used in your lab.
5. This is the key step: write a viper_env.bash script that looks like this:
   
```
export CONDA_ROOT=/home/viper/local/miniconda3  
export CONDA_ENVS_PATH=/home/viper/local/miniconda3/envs  
export PATH=$CONDA_ENVS_PATH/viper/bin:$PATH  
unset PYTHONPATH  
```
   
   - CONDA_ROOT is where you installed miniconda for user 'viper'  
   - CONDA_ENVS_PATH should simply be $CONDA_ROOT/envs  
   - PATH is overriding the user's PATH variable  
   - and the last command is to UNSET the user's PYTHONPATH, just in case they set it b/c we want them to use 'viper's CONDA PYTHONPATH.
6. So viper's home directory might look like this:  
   >config.yaml  
   >metasheet.csv  
   >viper/  
   >viper_env.bash
   
7. The final step is to write a simple bash script: viperSetup.sh to copy /home/viper/* to the local directory, i.e.:
   ```
   #!/bin/bash  
   cp -r /home/viper/* .
   ```
8. Finally, as sudo, place viperSetup.sh into a place like /usr/local/bin  

##Setting up and Running VIPER with the initialization procedure
So once you have setup a central viper user, the other users in the lab can 
start using viper by doing the following:

1. $ `mkdir PROJECT` #some arbitrary PROJECT name
2. $ `cd PROJECT`
3. PROJECT$ `viperSetup.sh` #which will copy the central viper to the local dir
4. PROJECT$ `source viper_env.bash` #which will take on the central viper's PATH
5. PROJECT$ `source activate viper` #which will activate the viper conda env  
\#This is where the user will have to define config.yaml and metasheet.csv  
6. (viper) PROJECT$ `snakemake -s viper/viper.snakefile -n` #to invoke a DRY-RUN
7. (viper) PROJECT$ `snakemake -s viper/viper.snakefile` #to invoke VIPER

##APPENDIX D: Running VIPER in a cluster
VIPER can be run on cluster systems which use a common filesystem for all compute-nodes (e.g. Sun Grid Engine) using the native [snakemake support](http://snakemake.readthedocs.io/en/latest/executable.html#cluster-execution) for these systems.  A separate "--cluster" parameter is given to the snakemake call to pass-along all of the relevant qsub commands, e.g.:
`snakemake -s viper/viper.snakefile -cluster "qsub {threads}" -j 32`

where -j determines the number of jobs to submit.

Please see your local cluster documentation as well as the snakemake documentation for more information.

##APPENDIX X: installing and running CDR3 analysis with TRUST
LEN TODO:
