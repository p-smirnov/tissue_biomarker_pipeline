library(data.table)

method <- "perm"



home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "smallPSets")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_perm_sig"))
myPvalDir <- file.path(".", paste0(method, "_meta_perm_sig"))


## First we load in the fixed effect permutation test results:


toRunPerm <- fread("toRunFiles/toRunMetaByGenePerm.txt", header=FALSE)

colnames(toRunPerm) <- c("Drug", 'Tissue', "Gene")

# toRunPerm[,hetFilePath := file.path(scratch, "perm_meta_hetPermTest", 
# 		paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_hetPerm_10000_out.rds"))]

# toRunPerm[,hetTestExists := file.exists(hetFilePath)]

# toRunPerm[hetTestExists == TRUE]


# toRunPerm[,PermFilePattern := paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_perm*")]

# toRunPerm[,PermFilePath := sapply(PermFilePattern, function(x) list.files(path=file.path(scratch, "perm_meta_bootstrap"), pattern=x, full.names=TRUE))]

# toRunPerm[,PermTestExists := sapply(PermFilePath, file.exists)]
# toRunPerm[!sapply(PermTestExists, any)]

toRunPerm[,PermResFile := file.path(myPvalDir, paste0("permSig_", make.names(Drug), "_", make.names(Tissue), "_", Gene,"_out.rds"))]
toRunPerm[,all(sapply(PermResFile, file.exists))]

system.time(toRunPerm[,PermRes := sapply(PermResFile, function(x) list(readRDS(x)))])

toRunPerm[,pvalue := sapply(PermRes, function(x) return(x[[1]]))]
toRunPerm[,estimate := sapply(PermRes, function(x) return(x[[2]]$t0))]

toRunPerm[, lower := sapply(PermRes, function(x) return(x[["ci"]]$bca[4]))]
toRunPerm[, upper := sapply(PermRes, function(x) return(x[["ci"]]$bca[5]))]


toRunPerm[order(pvalue)]

toRunPermOut <- toRunPerm

toRunPermOut <- toRunPermOut[,PermRes:=NULL]
toRunPermOut <- toRunPermOut[,PermResFile:=NULL]


write.csv(toRunPermOut, file="rnaResults/biomarker_res/BiomarkerPermRes.csv")


### Now lets do the bootstrap results



home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "smallPSets")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
bootSigDir <- file.path(scratch, "perm_meta_boot_sig")
bootSigDir <- file.path(".", paste0(method, "_meta_boot_sig"))



toRunBoot <- fread("toRunFiles/toRunMetaByGeneBoot.txt", header=FALSE)

colnames(toRunBoot) <- c("Drug", 'Tissue', "Gene")




toRunBoot[,BootSigPath := file.path(bootSigDir, paste0("bootSig_", make.names(Drug), "_", make.names(Tissue), "_", Gene, "_out.rds"))]

toRunBoot[,BootSigExists := sapply(BootSigPath, file.exists)]
toRunBoot[BootSigExists == FALSE]
toRunBoot <- toRunBoot[BootSigExists == TRUE]
## Job 1318428 is checking what is going on with these

toRunBoot <- toRunBoot[,BootRes := sapply(BootSigPath, function(x) return(list(readRDS(x))))]


toRunBoot[,pvalue := sapply(BootRes, function(x) return(x[[1]]))]
toRunBoot[,estimate := sapply(BootRes, function(x) return(x[[2]]))]

toRunBoot[, lower := sapply(BootRes, function(x) return(x[["ci"]]$bca[4]))]
toRunBoot[, upper := sapply(BootRes, function(x) return(x[["ci"]]$bca[5]))]


toRunBoot[order(pvalue)]
toRunBootOut <- toRunBoot

toRunBootOut <- toRunBootOut[,BootRes:=NULL]
toRunBootOut <- toRunBootOut[,BootSigExists:=NULL]
toRunBootOut <- toRunBootOut[,BootSigPath:=NULL]


write.csv(toRunBootOut, file="biomarker_res/BiomarkerBootRes.csv")

