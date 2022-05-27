# README


This repository contains scripts to run a meta-analysis of biomarkers from in vitro PharmacoGx datasets on HPC or the cloud. As all data is available, not just summary statistics for this analysis, we use mixed-linear models to come up with consensus estimates for the correlation between a gene expression and drug response across the in-vitro data, either pan-cancer of within a particular tissue type.  

**Unfortunately, the HPC pipeline does not currently work using the latest version of Snakemake, due to a bug in group-job dispatch. However, there are example scripts for running these jobs on Azure cloud services in the directory**

Due to some quirks on how Snakemake pulls and pushes files to Kubernetes working pods for running in the cloud, I have found it easier to keep all files in the same folder. 


## Preparing the pipeline

In general, the pipeline is run in three steps:

1. First, PharmacoSets are downloaded from orcestra.ca , and prepared for the pipeline. This is inevitably dataset specific, as choices
need to be made regarding which molecular type to use as well as how and whether to filter drug sensitivity data, and whether any other transformations
on the data should be made. An example of this step is found in the file: `makePSetsForRuns.R`. 
2. Then a file is created listing all the Gene-Tissue-Drug-Dataset combinations that will be evaluated by the pipeline. This is done so that throughout the 
run, there is a record of what was expected to be created, and what is actually output, helping catch any errors in the scripts or in execution (for example,
failed jobs). This is also the step where any filtering is done to minimum samples sizes per condition (20), or minimum number of responses. Any tuples that do not 
have sufficient data are excluded from this file. This script is expected to output a "masterToRunFile", with the structure as in "cnvCGSToRunList.txt". An 
example of this file is `makeMasterToRunRNA.R` . Note, that if pancancer analysis is desired, the pipeline takes `all` in the tissue column to indicate no 
filtering to tissue types. 
3. Finally, the pipeline can be run. The `Snakefile` runs the pipeline, each step of which is described in more detail below. It is configured using a YAML
file as in `CNV_cgc_cloud.yaml`. The configuration is also described below.


## Pipeline configuration 

```
master_table: # This is the name/path to the "masterToRunFile"
scratch_dir: # This is the path to the scratch directory, could be the same as the project directory. This is used to store 
             # intermediate files to pass data into and from the Julia scripts. Unfortunately, calling R from Julia does not 
             # work on my HPC environment. 
project_dir: # The name of the project directory 
data_dir: # The location of all the data. The psets are expected in data_dir/molDataType
molDataType: # The name for the molecular data type under examination. Used to create subfolders to organize results. 
max_first_stage_perm: # an optional cutoff for the maximum number of permutations to use in the first step
code_dir: # path to the code directory. On the cloud, this almost certainly should be ".", since snakemake pulls the 
          # code and works in the same directory. 
julia_working_dir: # The path to a directory which can function as $JULIAHOME. On some HPC, the unix home directory is not
                   # writable from compute nodes, so this needs to be modified to a writable directory. Also used as the 
                   # project directory for julia, to allow installing and precompiling packages for the script ahead of time. 

```

## Pipeline Structure 

The pipeline itself runs three main steps. 

1. First, the significance of each association (Drug, Tissue, Gene) triplet is evaluated in each dataset. For this, Snakemake reads the 
"masterToRunFile", and determines which Dataset-Tissue-Drug triplets exist. Then, it dispatches a job running "runDatasetSpecificPermutations.R" for 
that combination. This script will check which genes need to be evalated for this Dataset-Tissue-Drug combination, conduct an adaptive permutation test 
for each gene-drug association using the QUICK-STOP algorithm, and output a file with the results in a `PharmacoSig` object from `PharmacoGx`. Finally, 
once all jobs from this step are done, it runs the `makeToRunMetaByGene.R` file, which reads in the results from each of these jobs, and outputs a file that 
contains the significance call for each Drug-Tissue-Gene-Dataset tuple. This is then used by snakemake to determine which jobs should be dispatched in the next step.
2. The next step tests for heterogeneity in the effect between datasets, to determine if a fixed or random effect model should be used in the meta-analysis. For each 
Drug-Tissue-Gene triplet that had at least one significant association in one datasets from the previous step, it dispatches `runInterStudyHetPerm.R`, which writes 
out the results of a permutation test for inter-study heterogeneity to an `rds` file. Finally, it runs `makeToRunMetaByGeneForBoot.R`, which reads in all the results, 
and prints out a large text file collating the results of the heterogeneity test for each triplet. This again is used by Snakemake to decide which jobs to dispatch in 
the next step. 
3. Finally, the meta-analysis step can be run. Snakemake reads in the results of the previous step, and dispatches a job down one of two branches for each Gene-Drug-Tissue triplet, the fixed effect or random effect branch:
- Fixed effect branch: This is done in one step, Snakemake runs `runFixedEffectPerm.R` for each triplet, is an R script wrapping the permutations done in C code, either in `metaPermC.c`, or `metaPermCTissue.c`. Permutation testing is used to evaluate significance, and a bootstrap is then run to calculate 95% confidence intervals, primarily for plotting. 
- Random effect branch: This is more complicated. As mentioned above, on some HPC environments, it is tricky to set up R to call Julia (I suspect this has to do with Julia running in a singularity container). Therefore, first this branch writes out the data across all datasets for each Drug-Tissue-Gene triplet to a text file, either one by one `prepareBootDataForJulia.R` or in a batch `prepareBootDataForJuliaBatch.R`. It then calls `runMetaBootThreaded.jl` on each input file, which conducts a large bootstrap (1e7) samples for the purposes of bootstrap hypothesis testing. This script writes out results to a text file, which is read in by `evaluateJuliaResults.R` to calculate p values and confidence intervals for each meta-effect estimated. 

    The output of these steps ends with a file for Drug-Tissue-Gene triplet with the estimated average effect size across studies (the fixed effect from either the fixed or random effect models), the estimated p-value, and confidence interval bands.


## Running in the cloud

First, I recommend following the tutorial from the Snakemake authors for executing on kubernetes here: https://snakemake.readthedocs.io/en/stable/executing/cloud.html#generic-cloud-support-via-kubernetes

The files `cluster_setup.sh` and `runOnAKS.sh` are helper files for setting up specifically for Azure Kubernetes Service. 

## Running on HPC. 

A lot will depend on the configuration of your cluster. An example script for running jobs on the HPC4Health Slurm cluster at UHN is: `runOnH4H.sh`.

## Other scripts in repository

There are other helper scripts in the repository, such as `createMetaPlot*` for outputting forest plots of the meta-analysis for each marker (this is deprecated as a dashboard is under development to explore these results), and `assembleMetaResults.R` prints out tables with the results from the analysis collated into one file. 

