#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#----------------------------------------
# @authors: Tosh, Mahesh Vangala
# @emails: , vangalamaheshh@gmail.com
# @date: July, 1st, 2016
# @modified by Yuting Liu
# @modified date: Jun, 24, 2020
#----------------------------------------

from scripts.utils import _getProcessedCuffCounts

rule pca_plot:
    input:
        rpkmFile = "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.filtered_top{numgenes_plots}VarGenes.csv",
        annotFile = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
      #  expand("analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/images/pca_plot_{metacol}.png", metacol=config["metacols"]),
        pca_plot_out="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/pca_plot.pdf"
    params:
        pca_out_dir = "analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/"
    message: "Generating PCA plots"
    benchmark:
        "benchmarks/" + config["token"] + "/pca_plot_top{numgenes_plots}VarGenes.txt"
    shell:
        "mkdir -p {params.pca_out_dir}/images && source R-3.6.1.sh && Rscript ./modules/scripts/pca_plot.R {input.rpkmFile} {input.annotFile} {params.pca_out_dir}"

rule heatmapSS_plot:
    input:
        rpkmFile = "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.filtered_top{numgenes_plots}VarGenes.csv",
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        ss_plot_out="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSS_Spearman_plot.pdf",
        ss_txt_out="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSS_Spearman.txt",
        ss_plot_out2="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSS_Pearson_plot.pdf",
        ss_txt_out2="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSS_Pearson.txt"
    message: "Generating Sample-Sample Heatmap"
    params:
        ss_out_dir = "analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/",
        token=config['token'],
    benchmark:
        "benchmarks/" + config["token"] + "/heatmapSS_plot_top{numgenes_plots}VarGenes.txt"
    shell:
        "mkdir -p analysis/{params.token}/plots/images &&  source R-3.6.1.sh && Rscript ./modules/scripts/heatmapSS_plot.R {input.rpkmFile} "
        "{input.annotFile} {params.ss_out_dir} "


rule heatmapSF_plot:
    input:
        rpkmFile = "analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.filtered_top{numgenes_plots}VarGenes.csv",
        annotFile=config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        sf_plot_out="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSF_Spearman_plot.pdf",
        sf_txt_out="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSF_Spearman.txt",
        sf_plot_out2="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSF_Pearson_plot.pdf",
        sf_txt_out2="analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/heatmapSF_Pearson.txt"
    params:
        num_kmeans_clust = config["num_kmeans_clust"],
        sf_out_dir = "analysis/" + config["token"] + "/plots/top{numgenes_plots}VarGenes/",
        token=config['token'],
    message: "Generating Sample-Feature heatmap"
    benchmark:
        "benchmarks/" + config["token"] + "/heatmapSF_plot_top{numgenes_plots}VarGenes.txt"
    shell:
        "mkdir -p analysis/{params.token}/plots/images &&  source R-3.6.1.sh && Rscript ./modules/scripts/heatmapSF_plot.R {input.rpkmFile} "
        "{input.annotFile} {params.num_kmeans_clust} {params.sf_out_dir} "
