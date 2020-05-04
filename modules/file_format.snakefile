#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#------------------------------------


rule bam_to_bigwig:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        bai="analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    output:
        protected("analysis/bam2bw/{sample}/{sample}.CPM.bw")
    threads: 20
    message: "Converting {wildcards.sample} bam to bigwig"
    benchmark:
        "benchmarks/{sample}/{sample}.bam_to_bigwig.txt"
    shell:
        "bamCoverage --bam {input.bam} -o {output} "
        "--binSize 1 --normalizeUsing CPM --extendReads 0 -p {threads}"
