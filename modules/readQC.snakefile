#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#----------------------------------

rule read_distrib_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/read_distrib/{sample}.txt")
    message: "Running RseQC read distribution on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.read_distrib_qc.txt"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        python2=config['python2'],
        rseqc_path=config['rseqc_path'],
        bed_file=config['bed_file'],
    shell:
        "{params.pypath} {params.python2} {params.rseqc_path}/read_distribution.py"
        " --input-file={input}"
        " --refgene={params.bed_file} 1>{output}"

rule read_distrib_qc_matrix:
    input:
        read_distrib_files=expand( "analysis/RSeQC/read_distrib/{sample}.txt", 
            sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        matrix="analysis/" + config["token"] + "/RSeQC/read_distrib/read_distrib.matrix.tab",
        png="analysis/" + config["token"] + "/RSeQC/read_distrib/read_distrib.png"
    message: "Creating RseQC read distribution matrix"
    benchmark:
        "benchmarks/" + config["token"] + "/read_distrib_qc_matrix.txt"
    run:
        file_list_with_flag = " -f ".join( input.read_distrib_files )
        shell( "perl ./modules/scripts/read_distrib_matrix.pl -f {file_list_with_flag} 1>{output.matrix}" )
        shell( " source R-3.6.1.sh && Rscript ./modules/scripts/read_distrib.R {output.matrix} {output.png}" )


rule gene_body_cvg_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png"),
        protected("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r")
    threads: 8
    message: "Creating gene body coverage curves"
    benchmark:
        "benchmarks/{sample}/{sample}.gene_body_cvg_qc.txt"
    params: 
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        ds_bam = lambda wildcards: "analysis/RSeQC/gene_body_cvg/" + wildcards.sample + "/" + wildcards.sample + ".ds.bam",
        python2=config['python2'],
        rseqc_path=config['rseqc_path'],
        bed_file=config['bed_file'],
    shell:
        "picard DownsampleSam VALIDATION_STRINGENCY=LENIENT I={input}"
        " O={params.ds_bam} P=$(samtools flagstat -@ {threads} {input} |"
        " perl -e 'my $line = <STDIN>; chomp $line; my( $one, $two ) = ($line =~ /(\d+)\s+\+\s+(\d+)/); my $total = $one + $two;  my $one_M = $total < 1000000 ? 1 : (1000000 / $total); my $final_val = sprintf(\"%.2f\",$one_M); print $final_val > 0.00 ? $final_val : 0.01;') && "
        "samtools index {params.ds_bam} && "
        "{params.pypath} {params.python2} {params.rseqc_path}/geneBody_coverage.py -i {params.ds_bam} -r {params.bed_file}"
        " -f png -o analysis/RSeQC/gene_body_cvg/{wildcards.sample}/{wildcards.sample}"

rule plot_gene_body_cvg:
    input:
        samples_list=expand("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", 
            sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        rscript="analysis/" + config["token"] + "/RSeQC/gene_body_cvg/geneBodyCoverage.r",
        png="analysis/" + config["token"] + "/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        png_curves="analysis/" + config["token"] + "/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png"
    message: "Plotting gene body coverage"
    benchmark:
        "benchmarks/" + config["token"] + "/plot_gene_body_cvg.txt"
    shell:
        "perl ./modules/scripts/plot_gene_body_cvg.pl --rfile {output.rscript} --png {output.png} --curves_png {output.png_curves}"
        " {input.samples_list} &&  source R-3.6.1.sh && Rscript {output.rscript}"

rule junction_saturation:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf")
    message: "Determining junction saturation for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.junction_saturation.txt"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        python2=config['python2'],
        rseqc_path=config['rseqc_path'],
        bed_file=config['bed_file'],
    shell:
        "{params.pypath} {params.python2} {params.rseqc_path}/junction_saturation.py -i {input} -r {params.bed_file}"
        " -o analysis/RSeQC/junction_saturation/{wildcards.sample}/{wildcards.sample}"


rule collect_insert_size:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    params:
        picard_path=config['picard_path'],
        ref_fasta=config['ref_fasta'],
    output:
        protected("analysis/RSeQC/insert_size/{sample}/{sample}.histogram.pdf")
    message: "Collecting insert size for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.collect_insert_size.txt"
    shell:
        "{params.picard_path} CollectInsertSizeMetrics"
        " H={output} I={input} O=analysis/RSeQC/insert_size/{wildcards.sample}/{wildcards.sample} R={params.ref_fasta}"



