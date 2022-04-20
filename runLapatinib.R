## This code runs the first step of the pipeline, where adaptive permutation testing is used
## to determine if any gene is significantly correlated with DR within its own dataset.
## This is done with the help of functions from PharmacoGx 


library(PharmacoGx)
library(RhpcBLASctl)
library(data.table)
## making sure that if R was compiled to run on multiple cores, each spawned thread only uses 1 of them
print(RhpcBLASctl::blas_get_num_procs())
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
args <- commandArgs(trailingOnly = TRUE)

psetName <- args[[1]]

drug <- args[[2]]

tissue <- args[[3]]
print(paste("PSet:", psetName, "Drug:", drug, "Tissue:", tissue)) 
message(paste("PSet:", psetName, "Drug:", drug, "Tissue:", tissue)) 

nthread <- as.numeric(args[[4]])
myToRunFileName <- args[[5]]

if(!args[[6]]=="Default"){
	options("PharmacoGx_Max_Perm"=as.numeric(args[[6]]))
}



## Reading in paths from env variables

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")

dataDir <- Sys.getenv("DATA")

myDataDir <- dataDir
myOutDir <- file.path(project, "pearson_perm_res")
myRunDir <- file.path(project, "runlist_files")

containername <- Sys.getenv("containername", unset=NA_character_)




if(!is.na(containername)){
	myOutDir <- file.path(containername, myOutDir)
	myRunDir <- file.path(containername, myRunDir)
	myDataDir <- file.path(containername, myDataDir)
	project <- file.path(containername, project)
	scratch <- file.path(containername, scratch)


}



badchars <- "[ ]|[/]|[:]|[-]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))


toRun <- fread(myToRunFileName, header=FALSE)

colnames(toRun) <- c("Gene", "Tissue", "Drug", "PSet")

print(drug)
print(tissue)

# need to do this "trick" because names are made path safe, and arguments are derived from paths for snakemake's sake 
toRunThis <- toRun[make.names.2(toRun[,Drug]) == drug & make.names.2(toRun[,Tissue]) == tissue, ]
drug <- unique(toRunThis[,Drug])
tissue <- unique(toRunThis[,Tissue])
# pSets <- toRunThis[,4]

print(drug)
print(tissue)


## Loading in the dataset
switch(psetName, 
	   CCLE = {
	   		pset <- readRDS(file.path(myDataDir,"CCLE.rds"))
	   }, CCLE.CTRPv2 = {
	   		pset <- readRDS(file.path(myDataDir,"CCLE.CTRPv2.rds"))
	   }, CCLE.PRISM = {
	   		pset <- readRDS(file.path(myDataDir,"CCLE.PRISM.rds"))
	   }, GDSC_v1 = {
	   		pset <- readRDS(file.path(myDataDir,"GDSC1.rds"))
	   }, GDSC_v2 = {
	   		pset <- readRDS(file.path(myDataDir,"GDSC2.rds"))
	   }, gCSI = {
	   		pset <- readRDS(file.path(myDataDir,"gCSI.rds"))
	   }, GRAY = {
	   		pset <- readRDS(file.path(myDataDir,"GRAY.rds"))
	   }, UHNBreast = {
	   		pset <- readRDS(file.path(myDataDir,"UHNBreast.rds"))
	   }, Tavor = {
	   		pset <- readRDS(file.path(myDataDir, "Tavor.rds"))
	   }, BeatAML = {
	   		pset <- readRDS(file.path(myDataDir, "BeatAML.rds"))
	   }, "FIMM-AML-MCM" = {
	   		pset <- readRDS(file.path(myDataDir, "FIMM_MCM.rds"))
	   }, {stop("Please Provide a valid pset")})

## datasets used are subsetted to a single data type for efficiency. Maybe this should be passed in from config?
mData <- mDataNames(pset) 

## microarray and rnaseq annotations have different column names
gene_type_col <- ifelse("GeneBioType" %in% colnames(featureInfo(pset, mData)), "GeneBioType", "gene_type") 
## limiting feature space for power
ft <- rownames(featureInfo(pset, mData))[featureInfo(pset, mData)[[gene_type_col]] == "protein_coding"]



if(is.na(tissue) || tissue == "all"){
	chosen.cells <- cellNames(pset)
} else {
	chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
}
if(!length(chosen.cells)){
	stop("Something seems to have gone wrong with the provided tissue")
}

if(is.na(drug)){
	drug <- drugNames(pset)
}


# myToRunFileName <- file.path(myRunDir,"geneExpressionMasterToRunList.txt")

filteredFeatureList <- toRun[PSet == psetName, unique(Gene)]

ft <- ft[gsub(x=ft, pat="\\.[0-9]+$", rep="") %in% filteredFeatureList]

## run the permutation test for each gene in ft, for the drug and tissue selected. 
signature <- drugSensitivitySig(pset, mData, drugs=drug, features=ft,
	sensitivity.measure = "aac_recomputed", modeling.method="pearson", 
	inference.method="resampling", cells=chosen.cells, nthread=nthread, parallel.on = "gene")

if(!file.exists(myOutDir)){
	dir.create(myOutDir)
}

saveRDS(signature, file = file.path(myOutDir, make.names(paste0("signature_", psetName, "_", drug, "_", tissue, ".rds"))))
# ENSG00000000003,Bowel,5-Fluorouracil,CCLE.CTRPv2,

## after opt, without multicore correlations: 971.485 seconds

##CCLE.CTRPv2_PLX4720_Lymphoid.rds