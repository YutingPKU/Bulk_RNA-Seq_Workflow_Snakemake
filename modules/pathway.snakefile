#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-----------------------------------
# @authors: Tosh, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-----------------------------------

rule goterm_analysis:
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv",
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        out_file = "analysis/" + config["token"] + "/GO/{comparison}/{comparison}.DEGFDR{adjpvalcutoff}.goterm.done"
    params:
        gotermadjpvalcutoff = "{adjpvalcutoff}",
        numgoterms = config["numgoterms"],
        reference = config["reference"],
        prefix = "analysis/" + config["token"] + "/GO/{comparison}/{comparison}"
    message: "Creating Goterm Analysis plots for Differential Expressions for {wildcards.comparison}"
    benchmark:
        "benchmarks/" + config["token"] + "/{comparison}.DEGFDR{adjpvalcutoff}.goterm_analysis.txt"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/goterm_analysis.R {input.deseq} {params.gotermadjpvalcutoff} "
        "{params.numgoterms} {params.reference}  {params.prefix} && "
        " touch {output.out_file} "




