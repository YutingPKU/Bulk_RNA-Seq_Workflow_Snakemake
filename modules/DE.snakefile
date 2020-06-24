#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from scripts.utils import _getSTARcountsRes

rule limma_and_deseq:
    input:
        counts = _getSTARcountsRes(config)[0]
    output:
        limma = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.limma.csv",
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        deseqSum = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.sum.csv",
        #annotations
        limma_annot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.limma.annot.csv",
        deseq_annot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.annot.csv",

    params:
        s1 = lambda wildcards: ",".join(config['comps'][wildcards.comparison]['control']),
        s2 = lambda wildcards: ",".join(config['comps'][wildcards.comparison]['treat']),
        gene_annotation = config['gene_annotation'],
        prefix = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}"
    message: "Running differential expression analysis using limma and deseq for {wildcards.comparison}"
    benchmark:
        "benchmarks/" + config["token"] + "/{comparison}.limma_and_deseq.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/DEseq.R \"{input.counts}\" \"{params.s1}\" \"{params.s2}\" " 
        "{output.limma} {output.deseq} {output.limma_annot} {output.deseq_annot} "
        "{output.deseqSum} {params.gene_annotation} {params.prefix}"
        # ridiculous hack for singletons (Mahesh Vangala)
        "&& touch {output.limma} "
        "&& touch {output.limma_annot}"


rule deseq_limma_fc_plot:
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        limma = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.limma.csv"
    output:
        out_csv = "analysis/" + config["token"] + "/diffexp/{comparison}/deseq_limma_fc_corr.csv",
        out_png = "analysis/" + config["token"] + "/diffexp/{comparison}/deseq_limma_fc_corr.png"
    message: "Creatting deseq-limma correlation plot for {wildcards.comparison}"
    benchmark:
        "benchmarks/" + config["token"] + "/{comparison}.deseq_limma_fc_plot.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/deseq_limma_fc_corr.R {input.deseq} {input.limma} {output.out_csv} {output.out_png}"


rule fetch_DE_gene_list:
    input:
        deseq_file_list=expand("analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",comparison=config["comparisons"]),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/" + config["token"] + "/diffexp/de_summary.csv",
        png="analysis/" + config["token"] + "/diffexp/de_summary.png"
    message: "Creating Differential Expression summary"
    benchmark:
        "benchmarks/" + config["token"] + "/fetch_DE_gene_list.txt"
    run:
        deseq_file_string = ' -f '.join(input.deseq_file_list)
        shell("perl ./modules/scripts/get_de_summary_table.pl -f {deseq_file_string} 1>{output.csv}")
        shell(" source R-3.6.1.sh && Rscript ./modules/scripts/de_summary.R {output.csv} {output.png}")


#Generate volcano plots for each comparison
rule volcano_plot:
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        plot = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}_volcano.pdf",
        png = "analysis/" + config["token"] + "/plots/images/{comparison}_volcano.png"
    message: "Creating volcano plots for Differential Expressions for {wildcards.comparison}"
    benchmark:
        "benchmarks/" + config["token"] + "/{comparison}.volcano_plot.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/volcano_plot.R {input.deseq} {output.plot} {output.png}"


