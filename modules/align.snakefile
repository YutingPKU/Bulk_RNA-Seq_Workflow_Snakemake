#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab

#-------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

strand_command=""
rRNA_strand_command=""

if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
    rRNA_strand_command="--outFilterIntronMotifs RemoveNoncanonical"
else:
    rRNA_strand_command="--outSAMstrandField intronMotif"

run_fusion= True if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else False
gz_command="--readFilesCommand zcat" if config["samples"][config["ordered_sample_list"][0]][0][-3:] == '.gz' else ""


rule run_STAR:
    input:
        getFastq
    output:
        bam=protected("analysis/STAR/{sample}/{sample}.sorted.bam"),
        counts="analysis/STAR/{sample}/{sample}.counts.tab",
        log_file="analysis/STAR/{sample}/{sample}.Log.final.out",
    params:
        stranded=strand_command,
        gz_support=gz_command,
        prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        star_index = config['star_index'],
    threads: 20
    message: "Running STAR Alignment on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_STAR.txt"
    shell:
        "STAR --runMode alignReads --runThreadN {threads}"
        " --genomeDir {params.star_index}"
        " --readFilesIn {input} {params.gz_support}"
        " --outFileNamePrefix {params.prefix}."
        " --outSAMstrandField intronMotif"
        " --outSAMmode Full --outSAMattributes All {params.stranded}"
        " --outSAMattrRGline {params.readgroup}"
        " --outSAMtype BAM SortedByCoordinate"
        " --limitBAMsortRAM 50000000000" #50 Gib
        " --quantMode GeneCounts"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"

rule index_bam:
    """INDEX the {sample}.sorted.bam file"""
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    message: "Indexing {wildcards.sample}.sorted.bam"
    benchmark:
        "benchmarks/{sample}/{sample}.index_bam.txt"
    shell:
        "samtools index {input}"

rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
        star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/" + config["token"] + "/STAR/STAR_Align_Report.csv",
        png="analysis/" + config["token"] + "/STAR/STAR_Align_Report.png",
        gene_counts="analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv"
    message: "Generating STAR report"
    #priority: 3
    benchmark:
        "benchmarks/" + config["token"] + "/generate_STAR_report.txt"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl ./modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( " source R-3.6.1.sh && Rscript ./modules/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl ./modules/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule batch_effect_removal_star:
    input:
        starmat = "analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv",
        annotFile = config["metasheet"]
    output:
        starcsvoutput="analysis/" + config["token"] + "/STAR/batch_corrected_STAR_Gene_Counts.csv",
        starpdfoutput="analysis/" + config["token"] + "/STAR/star_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "star",
        token = config['token'],
    message: "Removing batch effect from STAR Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    #priority: 2
    benchmark:
        "benchmarks/" + config["token"] + "/batch_effect_removal_star.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/batch_effect_removal.R {input.starmat} {input.annotFile} "
        "{params.batch_column} {params.datatype} {output.starcsvoutput} {output.starpdfoutput} "
        " && mv {input.starmat} analysis/{params.token}/STAR/without_batch_correction_STAR_Gene_Counts.csv"



rule align_SJtab2JunctionsBed:
    """Convert STAR's SJ.out.tab to (tophat) junctions.bed BED12 format"""
    input:
        "analysis/STAR/{sample}/{sample}.SJ.out.tab"
    output:
        "analysis/STAR/{sample}/{sample}.junctions.bed"
    benchmark:
        "benchmarks/{sample}/{sample}.align_SJtab2JunctionsBed.txt"
    shell:
        "./modules/scripts/STAR_SJtab2JunctionsBed.py -f {input} > {output}"
