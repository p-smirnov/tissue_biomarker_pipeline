snakemake --cluster "sbatch -J {params.jobname} -o /scratch/b/bhaibeka/psmirnov/slurm-%j.log -t {params.runtime} -c 40 " \
 --groups hetTest=group1 evaluateJuliaRes=group4 --group-components group1=30 group4=40 -s NiagaraSnakefile\
  --jobs 1 --latency-wait 100 --scheduler greedy --rerun-incomplete all