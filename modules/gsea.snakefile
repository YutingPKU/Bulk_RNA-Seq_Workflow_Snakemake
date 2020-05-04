#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def gseaInputFn(wildcards):
    """We need to handle whether this is a mouse run--if so, convert
    the mouse genes to human genes in prep for GSEA analysis"""
    if 'assembly' in config and config['assembly'].startswith("mm"):
        return ["analysis/" + config["token"] + "/diffexp/%s/%s.deseq.converted.csv" % (wildcards.comparison, wildcards.comparison)]
    else:
        return ["analysis/" + config["token"] + "/diffexp/%s/%s.deseq.csv" % (wildcards.comparison, wildcards.comparison)]

rule convertMouseToHuman:
    """A rule to convert mouse to human"""
    input:
        deseq = "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.csv"
    params:
        map_f = "static/gsea/mmToHsGeneMap.csv"
    message: "GSEA: converting mouse genes to human genes"
    output:
        "analysis/" + config["token"] + "/diffexp/{comparison}/{comparison}.deseq.converted.csv"
    shell:
        "./modules/scripts/mmGeneToHs.py -m {params.map_f} -i {input} -o {output}"

rule gsea:
    """Use clusterProfile enchricher fn to perform gene set enrichment analysis
    using MSigDB (default) or user defined db. Uses pvalAdj cut-off of 0.05 to 
    determine gene list.
    """
    input:
        gseaInputFn
    output:
        gene_list = "analysis/" + config["token"] + "/gsea/{comparison}/{comparison}.gene_list.txt",
        gsea = "analysis/" + config["token"] + "/gsea/{comparison}/{comparison}.gene_set.enrichment.txt",
        dotplot = "analysis/" + config["token"] + "/gsea/{comparison}/{comparison}.gene_set.enrichment.dotplot.png",
    params:
        db = config["gsea_db"],
        out_path = "analysis/" + config["token"] + "/gsea/{comparison}/{comparison}",
        title = "{comparison}"
    message: "Running Gene Set Enrichment Analysis on {wildcards.comparison}"
    shell:
        " source R-3.6.1.sh && Rscript ./modules/scripts/gsea.R {input} {params.db} \"{params.title}\" {params.out_path}"


