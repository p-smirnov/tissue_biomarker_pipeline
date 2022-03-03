configfile: "CNV_cgc_cloud.yaml"


#TODO: set env variables through prefixes to commands
import pandas as pd
import os

print(os.getcwd())

masterToRunFile = config["master_table"]



def cUFL(fn):
    fn = [x.replace(' ', '.') for x in fn]
    fn = [x.replace('/', '.') for x in fn]
    fn = [x.replace(':', '.') for x in fn]
    fn = [x.replace('-', '.') for x in fn]
    fn = [x.replace('/', '.') for x in fn]
    fn = [x.replace(',', '.') for x in fn]
    fn = [x.replace('(', '.') for x in fn]
    fn = [x.replace(')', '.') for x in fn]
    fn = [x.replace('[', '.') for x in fn]
    fn = [x.replace(']', '.') for x in fn]
    fn = [x.replace('+', '.') for x in fn]
    fn = [x.replace("'", '.') for x in fn]
    return fn

def cUF(fn):
    fn = fn.replace(' ', '.')  
    fn = fn.replace('/', '.') 
    fn = fn.replace(':', '.') 
    fn = fn.replace('-', '.') 
    fn = fn.replace('/', '.')
    fn = fn.replace(',', '.')
    fn = fn.replace('(', '.')  
    fn = fn.replace(')', '.')    
    fn = fn.replace('[', '.')  
    fn = fn.replace(']', '.')
    fn = fn.replace('+', '.')
    fn = fn.replace("'", '.')
    return fn


scratch_dir = config["scratch_dir"]
project_dir = config["project_dir"]
data_dir = config["data_dir"]
mol_data_type = config["molDataType"]
code_dir = config["code_dir"]
julia_working_dir = config["julia_working_dir"]


scratch_dir = scratch_dir + "/" + mol_data_type
project_dir = project_dir + "/" + mol_data_type
data_dir = data_dir + "/" + mol_data_type

signature_dir = project_dir + "/pearson_perm_res"
runlist_dir = project_dir + "/runlist_files"
het_dir = project_dir + "/hetTest"
julia_data_dir = scratch_dir + "/data4juliaBoot"
julia_out_dir = scratch_dir + "/juliaBoot"
perm_out_dir = project_dir + "/perm_meta_out"
boot_out_sig = project_dir + "/perm_meta_boot_sig"
perm_out_sig = project_dir + "/perm_meta_perm_sig"


pset_list = ['Data/TBPInputs/cnv/CCLE.CTRPv2.rds', 'Data/TBPInputs/cnv/gCSI.rds', 'Data/TBPInputs/cnv/GDSC2.rds']
print(pset_list)

checkpoint getFirstStageToRun:
    input: masterToRunFile
    output: masterToRunFile + "2"
    params:
        runtime="1:00:00"
    resources:
        time="1:00:00"
    shell:
        """
        #! /bin/bash

        cp {input} {output}
        """
print(signature_dir)

def getFirstStageOutput(wildcards):
    with checkpoints.getFirstStageToRun.get(**wildcards).output[0].open() as f:
        master_table  = pd.read_csv(f, header=None)
        # pd.read_csv("/cluster/projects/bhklab/Projects/Tissue_Biomarker_Meta/runlist_files/geneExpressionMasterToRunList.txt", header=None)
        fst = master_table.drop_duplicates([1,2,3])[[1,2,3]]
        first_stage_ids = ['signature_{PSet}_{Drug}_{Tissue}.rds'.format(PSet = x[3], Tissue = x[1], Drug = x[2]) for (row,x) in fst.iterrows()]
        first_stage_ids = cUFL(first_stage_ids)
        first_stage_ids = [signature_dir + "/" + x for x in first_stage_ids]
        fst['ids'] = first_stage_ids
        fst = fst.set_index("ids")
        return(first_stage_ids)



checkpoint getSigGenes:
    input: getFirstStageOutput, toRunFile=rules.getFirstStageToRun.output
    output: runlist_dir + "/toRunMetaByGene.txt"
    params:
        runtime="24:00:00"
    resources: 
        time="24:00:00"
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0

        MKL_NUM_THREADS=1 MKL_DOMAIN_NUM_THREADS=1 OMP_NUM_THREADS=1\
         PROJECT={project_dir} DATA={data_dir} Rscript {code_dir}/makeToRunMetaByGene.R {input.toRunFile}
        """

rule runFirstStage:
    input: pset_list
    params:
        runtime="24:00:00", 
        # pset=lambda wildcards: fst.loc()[wildcards.id][3],
        # drug=lambda wildcards: fst.loc()[wildcards.id][2],
        # tissue=lambda wildcards: fst.loc()[wildcards.id][1],
        rthreads=20,
        max_number_of_perms=config["max_first_stage_perm"]
    threads: 1
    resources:
        time="24:00:00"
    output: signature_dir + "/signature_{PSet}_{Drug}_{Tissue}.rds"
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0
        SCRATCH={scratch_dir} 

        MKL_NUM_THREADS=1 MKL_DOMAIN_NUM_THREADS=1 OMP_NUM_THREADS=1 CODE={code_dir}\
         PROJECT={project_dir} DATA={data_dir} Rscript {code_dir}/runLapatinib.R {wildcards.PSet} "{wildcards.Drug}" '{wildcards.Tissue}' {threads} {masterToRunFile} {params.max_number_of_perms} 
        """

def outputFromHetTest(wildcards):
    with checkpoints.getSigGenes.get(**wildcards).output[0].open() as f:
        fstRes = pd.read_csv(f, header=None)
        fstRes = fstRes.loc()[fstRes[4]==1].drop_duplicates([0,1,2])[[0,1,2]]
        outFileNames = ['{output_dir}/{Drug}_{Tissue}_{Gene}_hetPerm_10000_out.rds'.format(Gene = cUF(x[0]), Tissue = cUF(x[1]), Drug = cUF(x[2]), output_dir=het_dir) for (row,x) in fstRes.iterrows()]
        # fstRes['filename'] = outFileNames
        # fstRes['ids'] = cUFL(['{Drug}_{Tissue}_{Gene}'.format(Gene = x[0], Tissue = x[1], Drug = x[2]) for (row,x) in fstRes.iterrows()])
        # fstRes = fstRes.set_index("ids")
        return outFileNames


rule hetTest:
    input: pset_list
    output: "{output_dir}/{Drug}_{Tissue}_{Gene}_hetPerm_10000_out.rds"
    params:
        R=10000
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0
        SCRATCH={scratch_dir} 


        MKL_NUM_THREADS=1 MKL_DOMAIN_NUM_THREADS=1 OMP_NUM_THREADS=1 PROJECT={project_dir} DATA={data_dir}\
         Rscript {code_dir}/runInterStudyHetPerm.R  {wildcards.Drug} {wildcards.Tissue} {wildcards.Gene} {params.R} > $SCRATCH/InterStudy.out
        """

checkpoint evaluteHeterogeneity:
    input: outputFromHetTest, runlist_dir + "/toRunMetaByGene.txt", pset_list
    output: runlist_dir + "/metaHetTestRes.txt"
    params:
        outdir = runlist_dir, 
        indir = runlist_dir
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0

        MKL_NUM_THREADS=1 MKL_DOMAIN_NUM_THREADS=1 OMP_NUM_THREADS=1 PROJECT={project_dir} DATA={data_dir}\
         Rscript {code_dir}/makeToRunMetaByGeneForBoot.R {params.indir} {params.outdir}
        """

checkpoint outputJuliaData:
    input:  runlist_dir + "/metaHetTestRes.txt", runlist_dir + "/toRunMetaByGene.txt", pset_list
    output: directory(julia_data_dir)
    params:
        runtime="24:00:00",
        outdir = julia_data_dir, 
        indir = runlist_dir
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0
        SCRATCH={scratch_dir} 


        MKL_NUM_THREADS=1 MKL_DOMAIN_NUM_THREADS=1 OMP_NUM_THREADS=1 PROJECT={project_dir} DATA={data_dir} \
        Rscript {code_dir}/prepareBootDataForJuliaBatch.R  {params.indir} {params.outdir}
        """


rule runJuliaH4H:
    input: julia_data_dir + "/modelData_{Gene}_{Drug}_{Tissue}.csv"
    output: "{output_dir}/metaBootRes_{Gene}_{Drug}_{Tissue}_out.txt"
    params:
        runtime="6:00:00",
        partition="all",
        R=10000, 
        juliathread= 1, 
    shell:
        """
        #! /bin/bash
        set +u;
        # module load NiaEnv/2019b julia/1.5.3
        SCRATCH={scratch_dir} 

        JULIA_DEPOT_PATH={julia_working_dir}/.julia julia --project={julia_working_dir} --threads {params.juliathread} \
        {code_dir}/runMetaBootThreaded.jl {wildcards.Drug} {wildcards.Tissue} {wildcards.Gene} {input} {wildcards.output_dir}
        """

rule evaluateJuliaRes:
    input: [julia_out_dir + "/metaBootRes_{Gene}_{Drug}_{Tissue}_out.txt", julia_data_dir + "/modelData_{Gene}_{Drug}_{Tissue}.csv"]
    output: "{output_dir}/bootSig_{Gene}_{Drug}_{Tissue}_out.rds"
    params: # todo fix these
        runtime="1:00:00",
        partition="all",
        R=10000,
        juliaDataDir=julia_data_dir
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0
        SCRATCH={scratch_dir} 

        PROJECT={project_dir} DATA={data_dir} \
        Rscript {code_dir}/evaluateJuliaResults.R {wildcards.Drug} {wildcards.Tissue} {wildcards.Gene} 10000 {julia_data_dir}/modelData_{wildcards.Gene}_{wildcards.Drug}_{wildcards.Tissue}.csv {params.juliaDataDir}/modelData_{wildcards.Gene}_{wildcards.Drug}_{wildcards.Tissue}.csv
        """


def requiredResJulia(wildcards):
    with checkpoints.evaluteHeterogeneity.get(**wildcards).output[0].open() as f:
        # test = checkpoints.outputJuliaData.get(**wildcards).output[0]
        # print(test)
        hetRes = pd.read_csv(f)
        hetRes = hetRes.loc()[hetRes["hetTestRes"]==1]
        outFileNames = ['{output_dir}/bootSig_{Gene}_{Drug}_{Tissue}_out.rds'.format(Gene = cUF(x[0]), Tissue = cUF(x[1]), Drug = cUF(x[2]), output_dir=boot_out_sig) for (row,x) in hetRes.iterrows()]
        return outFileNames


def requiredResFixedEffects(wildcards):
    with checkpoints.evaluteHeterogeneity.get(**wildcards).output[0].open() as f:
        # test = checkpoints.outputJuliaData.get(**wildcards).output[0]
        # print(test)
        hetRes = pd.read_csv(f)
        hetRes = hetRes.loc()[hetRes["hetTestRes"]==0]
        outFileNames = ['{output_dir}/permSig_{Gene}_{Drug}_{Tissue}_out.rds'.format(Gene = cUF(x[0]), Tissue = cUF(x[1]), Drug = cUF(x[2]), output_dir=perm_out_sig) for (row,x) in hetRes.iterrows()]
        return outFileNames


rule runFixedEffectTest:
    input: runlist_dir + "/metaHetTestRes.txt", runlist_dir + "/toRunMetaByGene.txt", pset_list
    output: [perm_out_sig + "/permSig_{Gene}_{Drug}_{Tissue}_out.rds", perm_out_dir + "/metaPermRes_{Gene}_{Drug}_{Tissue}_perm.RDS"]
    threads: 1
    params: #TODO fix me
        runtime="0:30:00",
        partition="all",
        rundir=runlist_dir,
        fixedeffthread=1,
        alpha=0.05
    resources:
        time="0:30:00", 
        mem_mb=7500
    shell:
        """
        #! /bin/bash
        set +u;
        # module load CCEnv  nixpkgs/16.09 gcc/8.3.0 r/4.0.0
        SCRATCH={scratch_dir} 
        R CMD SHLIB {code_dir}/metaPermC.c
        R CMD SHLIB {code_dir}/metaPermCTissue.c
       

        PROJECT={project_dir} DATA={data_dir} Rscript {code_dir}/runFixedEffectPerm.R \
        {wildcards.Drug} {wildcards.Tissue} {wildcards.Gene} {params.alpha} {threads} {params.rundir} {code_dir}
        """
def requiredResAll(wildcards):
    return requiredResJulia(wildcards) + requiredResFixedEffects(wildcards)


rule all:
        input: requiredResAll
        
