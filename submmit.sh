#!/bin/bash

#-------------------
# @author: Yuting Liu
# @email: lyt17@pku.edu.cn
# @date: Jun, 24, 2020
#--------------------
# how to run the workflow on cluster
 

#SLURM_ARGS="pkubatch -p cn-short -N 1 -c 20 --qos=lch3000cns -A lch3000_g1   -J {cluster.job-name} -o {cluster.output} -e {cluster.error}"
# print the running plan and command, not run 
snakemake -np -s viper.snakefile

# run the pipeline on cluster
nohup snakemake -j 10  -pr  -c "pkubatch -p cn-short -N 1 -c 20 --qos=lch3000cns -A lch3000_g1   -J {rule}.{wildcards} -o logs/cluster/{rule}/{rule}.{wildcards}_%j.out -e logs/cluster/{rule}/{rule}.{wildcards}_%j.err " -s viper.snakefile -k 2> snakemake.log &

# before rerun the pipeline, unlock the working directory
snakemake --unlock -np -s viper.snakefile
