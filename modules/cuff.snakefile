#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#--------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#--------------------------------

from scripts.utils import _getCuffCounts

cuff_command=""

if( config["stranded"] ):
    cuff_command="--library-type " + config["library_type"]


rule run_cufflinks:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        genes_cuff_out = protected("analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking"),
        iso_cuff_out = protected("analysis/cufflinks/{sample}/{sample}.isoforms.fpkm_tracking")
    threads: 20
    message: "Running Cufflinks on {wildcards.sample}"
    params:
        library_command=cuff_command,
        gtf_file=config['gtf_file'],
    benchmark:
        "benchmarks/{sample}/{sample}.run_cufflinks.txt"
    shell:
        "cufflinks -o analysis/cufflinks/{wildcards.sample} -p {threads} -G {params.gtf_file} {params.library_command} {input}"
        " && mv analysis/cufflinks/{wildcards.sample}/genes.fpkm_tracking {output.genes_cuff_out}"
        " && mv analysis/cufflinks/{wildcards.sample}/isoforms.fpkm_tracking {output.iso_cuff_out}"

rule generate_cuff_matrix:
    input:
        cuff_gene_fpkms=expand( "analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv"
    message: "Generating expression matrix using cufflinks counts"
    #priority: 1
    benchmark:
        "benchmarks/" + config["token"] + "/generate_cuff_matrix.txt"
    run:
        fpkm_files= " -f ".join( input.cuff_gene_fpkms )
        shell( "perl ./modules/scripts/raw_and_fpkm_count_matrix.pl -c -d -f {fpkm_files} 1>{output}" )


rule generate_gct_file:
    input:
        rpkmFile = _getCuffCounts(config)[1],
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.gct"
    message: "Outputting .GCT file from Raw Cuff Gene Counts File"
    #priority: 1
    benchmark:
        "benchmarks/" + config["token"] + "/generate_gct_file.txt"
    run:
        shell( " source R-3.6.1.sh && Rscript ./modules/scripts/GeneCountsToGCT.R {input.rpkmFile} {output}" )

rule generate_cuff_isoform_matrix:
    input:
        cuff_gene_fpkms=expand( "analysis/cufflinks/{sample}/{sample}.isoforms.fpkm_tracking", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        "analysis/" + config["token"] + "/cufflinks/Cuff_Isoform_Counts.csv"
    message: "Generating expression matrix using cufflinks isoform counts"
    #priority: 3
    params:
        #What to call our col 0
        iid="Transcript_ID"
    benchmark:
        "benchmarks/" + config["token"] + "/generate_cuff_isoform_matrix.txt"
    run:
        fpkm_files= " -f ".join(input.cuff_gene_fpkms)
        shell("./modules/scripts/cuff_collect_fpkm.py -n {params.iid} -f {fpkm_files} > {output}")

rule batch_effect_removal_cufflinks:
    input:
        cuffmat = "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv",
        annotFile = config["metasheet"]
    output:
        cuffcsvoutput="analysis/" + config["token"] + "/cufflinks/batch_corrected_Cuff_Gene_Counts.csv",
        cuffpdfoutput="analysis/" + config["token"] + "/cufflinks/cuff_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "cufflinks",
        token=config['token'],
    message: "Removing batch effect from Cufflinks Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    #priority: 2
    benchmark:
        "benchmarks/" + config["token"] + "/batch_effect_removal_cufflinks.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/batch_effect_removal.R {input.cuffmat} {input.annotFile} {params.batch_column} "
        "{params.datatype} {output.cuffcsvoutput} {output.cuffpdfoutput} "
        " && mv {input.cuffmat} analysis/{params.token}/cufflinks/without_batch_correction_Cuff_Gene_Counts.csv "


rule fpkm_plot:
    input:
        cuffmat = _getCuffCounts(config)[1],
        annotFile = config["metasheet"]
    output:
        fpkm_png = "analysis/" + config["token"] + "/plots/gene_counts.fpkm.png"
    message: "Plot gene counts at various fpkm cutoffs"
    benchmark:
        "benchmarks/" + config["token"] + "/fpkm_plot.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/fpkm_plot.R {input.cuffmat} {output.fpkm_png}"
