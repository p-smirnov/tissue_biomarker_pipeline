library(PharmacoGx)
library(data.table)
library(qs)

method <- "perm"

args <- commandArgs(trailingOnly = TRUE)
myToRunFileName <- args[[1]]

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")


dataDir <- Sys.getenv("DATA")


myDataDir <- dataDir

mySigDir <- file.path(project, paste0("pearson_",method,"_res"))
myRunDir <- file.path(project, "runlist_files")

## dealing with possibly running in the cloud
containername <- Sys.getenv("containername", unset=NA_character_)

if(!is.na(containername)){
	mySigDir <- file.path(containername, mySigDir)
	myRunDir <- file.path(containername, myRunDir)
}



# myToRunFileName <- file.path(myRunDir,"geneExpressionMasterToRunList.txt")

toRun <- fread(myToRunFileName, header=FALSE)

colnames(toRun) <- c("Gene", "Tissue", "Drug", "PSet")

toRunFirstStage <- unique(toRun[,.(PSet, Drug, Tissue)])

toRunFirstStage[,sig.file.path := file.path(mySigDir, make.names(paste0("signature_", PSet, "_", Drug, "_", Tissue, ".qs")))]

toRunFirstStage[,sig.obj := lapply(sig.file.path, qread)]

# toRunFirstStageRes <- toRunFirstStage[,.(Gene = rownames(sig.obj[[1]]), 
# 	Significant = sig.obj[[1]][,1,"significant"],
# Estimate = sig.obj[[1]][,1,"estimate"]), .(PSet, Drug, Tissue)]


toRunFirstStageRes <- toRunFirstStage[,.(Gene = rownames(sig.obj[[1]]), 
	Significant = sig.obj[[1]][,1,"significant"]), .(PSet, Drug, Tissue)]

toRunFirstStageRes[,Gene := gsub(Gene, pat="\\..+", rep="")]

toRunFirstStageRes <- merge(toRun, toRunFirstStageRes, by=c("Gene", "Tissue", "Drug", "PSet"), all.x=TRUE)

fwrite(toRunFirstStageRes, file=file.path(myRunDir,"toRunMetaByGene.txt"), sep=",", col.names=FALSE, quote=FALSE)

