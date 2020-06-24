#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
# @modified by Yuting Liu
# @modified date: Jun, 24, 2020
#-----------------------------------

def getTargetInfo(config):
    targetFiles = []
    targetFiles.extend([_getSTARaligns(config),
                        _getSTARcounts(config),
                        _getCuffCounts(config),
                        _getCuffIsoCounts(config),
                        _getProcessedCuffCounts(config),
                        _bw(config),
                        _DE(config),
                        _cluster(config),
                        _pathway(config),
                        _copyMetaFiles(config)])
    return targetFiles

def _getSTARaligns(config):
    """ensure that bam and indexes are built"""
    #STAR alignment sorted.bam file, its index, and transcript count
    ls = []
    if config["align"] == True:
        for sample in config["ordered_sample_list"]:
            ls.append("analysis/STAR/"+sample+"/"+sample+".sorted.bam")
            ls.append("analysis/STAR/"+sample+"/"+sample+".sorted.bam.bai")
        ls.append("analysis/" + config['token'] + "/STAR/STAR_Align_Report.png")
    return ls

## Returns proper count files for with and without batch effect correction
def _getSTARcounts(config):
    STAR_out_files = []
    if config["star_counts"] == True:
        STAR_out_files = ["analysis/" + config["token"] + "/STAR/batch_corrected_STAR_Gene_Counts.csv"] if config["batch_effect_removal"] == True else ["analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv"]
    return STAR_out_files

## Returns proper count files for with and without batch effect correction
def _getSTARcountsRes(config):
    STAR_out_files = ["analysis/" + config["token"] + "/STAR/batch_corrected_STAR_Gene_Counts.csv"] if config["batch_effect_removal"] == True else ["analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv"]
    return STAR_out_files


def _getCuffCounts(config):
    cuff_files = []
    if config["cuff_counts"] == True:
        cuff_files = ["analysis/" + config["token"] + "/plots/gene_counts.fpkm.png"]
        if config["batch_effect_removal"] == True:
            cuff_files.append("analysis/" + config["token"] + "/cufflinks/batch_corrected_Cuff_Gene_Counts.csv")
        else:
            cuff_files.append("analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.csv")
    return cuff_files

def _getCuffIsoCounts(config):
    cuff_files = []
    if config["cuff_counts"] == True:
        cuff_files = ["analysis/" + config["token"] + "/cufflinks/Cuff_Isoform_Counts.csv"]
    return cuff_files

def _getProcessedCuffCounts(config):
    ls = []
    if config["cuff_filter_counts"] == True:
      #  print(config["numgenes_plots"])
      #  print(type(config["numgenes_plots"]))
      #  nls = config["numgenes_plots"].replace("[","")
      #  nls = nls.replace("]","")
      #  nls = nls.replace(" ","")
        nls = config["numgenes_plots"].split(",")
      #  print(nls)
      #  print(type(nls))
        for numgenes_plots in nls:
            ls.append("analysis/" + config["token"] + "/cufflinks/Cuff_Gene_Counts.filtered_top" + str(numgenes_plots) +"VarGenes.csv" )
    return ls


def _bw(config):
    bw_files = []
    if config["bam2bw"]:
        bw_files.extend(["analysis/bam2bw/" + sample + "/" + sample + ".CPM.bw" for sample in config["ordered_sample_list"]])
    return bw_files



## Return the sample QC: PCA and clustering heatmap
def _cluster(config):
    cluster_files = []
    if config["cluster"] == True:
        nls = config["numgenes_plots"].split(",")
        for num in nls:
            num = str(num)
            cluster_files.append("analysis/" + config["token"] + "/plots/top"+ num + "VarGenes/pca_plot.pdf")
            cluster_files.append("analysis/" + config["token"] + "/plots/top" + num + "VarGenes/heatmapSS_Spearman_plot.pdf")
            cluster_files.append("analysis/" + config["token"] + "/plots/top" + num + "VarGenes/heatmapSF_Spearman_plot.pdf")
            cluster_files.append("analysis/" + config["token"] + "/plots/top" + num + "VarGenes/heatmapSS_Pearson_plot.pdf")
            cluster_files.append("analysis/" + config["token"] + "/plots/top" + num + "VarGenes/heatmapSF_Pearson_plot.pdf")
    return cluster_files

## Return the DEG using DESeq2 or edgeR
def _DE(config):
    de_list = []
    if config["comparisons"] and config["deg"] == True:
        de_list.append("analysis/" + config["token"] + "/diffexp/de_summary.png")
        de_list.extend([["analysis/" + config["token"] + "/diffexp/" + comp + "/" + comp + "_volcano.pdf",
                        "analysis/" + config["token"] + "/diffexp/" + comp + "/deseq_limma_fc_corr.png"]
            if len(config['comps'][comp]['control']) > 1 and len(config['comps'][comp]['treat']) > 1 else
            ["analysis/" + config["token"] + "/diffexp/" + comp + "/" + comp + "_volcano.pdf"] for comp in config["comparisons"]])
    return de_list


## Return the functional enrichment analysis for DEG
def _pathway(config):
    path_files = []
    if config["go_analysis"] == True:
        nls = config["goterm_adjpval_cutoff"].split(",")
        for comp in config["comparisons"]:
            for pval in nls:
                path_files.append("analysis/" + config["token"] + "/GO/" + comp + "/" + comp + ".DEGFDR" + pval + ".goterm.done")
    return path_files


def _copyMetaFiles(config):
    return ["analysis/" + config["token"] + "/" + config["token"] + '.config.yaml',
            "analysis/" + config["token"] + "/" + config["token"] + '.metasheet.csv']



