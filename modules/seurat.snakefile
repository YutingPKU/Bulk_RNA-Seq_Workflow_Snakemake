#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#-------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 25, 2016
#-------------------------

rule run_seurat:
    input:
        gene_matrix = "analysis/" + config["token"] + "/RSEM/tpm_gene_matrix.csv",
        metasheet = config["metasheet"]
    output:
        seurat_out = "analysis/" + config["token"] + "/seurat/seurat.done",
        seurat_out_dir = "analysis/" + config["token"] + "/seurat/"
    message: "Running Seurat - tSNE"
    benchmark:
        "benchmarks/" + config["token"] + "/run_seurat.txt"
    shell:
        "{config['seurat_path']}/Rscript --default-packages=stats,graphics,grDevices,utils,datasets,base,methods viper/modules/scripts/run_seurat.R "
        "{input.gene_matrix} {input.metasheet} {output.seurat_out_dir} "
        "&& touch {output.seurat_out}"
