SLURM_ARGS="pkubatch -p cn-short -N 1 -c 20 --qos=lch3000cns -A lch3000_g1   -J {cluster.job-name} -o {cluster.output} -e {cluster.error}"
snakemake -j 10  -pr --cluster-config cluster.yaml --cluster "$SLURM_ARGS" -s viper.snakefile  2> snakemake.log 
nohup snakemake -j 10  -pr  -c "pkubatch -p cn-short -N 1 -c 20 --qos=lch3000cns -A lch3000_g1   -J {rule}.{wildcards} -o logs/cluster/{rule}.{wildcards}_%j.out -e logs/cluster/{rule}.{wildcards}_%j.err " -s viper.snakefile -k 2> snakemake.log &
