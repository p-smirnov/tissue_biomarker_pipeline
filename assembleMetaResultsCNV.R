library(data.table)
corP <- 60671

## First we load in the fixed effect permutation test results:
myPvalDir <- "cnvRes/cnv/perm_meta_perm_sig"

hetTestRes <- fread("cnvRes/cnv/runlist_files/metaHetTestRes.txt", header=T)

permRes <- hetTestRes[hetTestRes==0]


bootSigDir <- "cnvRes/cnv/perm_meta_boot_sig"

toRunBoot <- hetTestRes[hetTestRes==1]



toRunBoot[,BootSigPath := file.path(bootSigDir, make.names(paste0("bootSig_", Gene,  "_", Drug, "_", Tissue,  "_out.rds")))]

toRunBoot[,BootSigExists := sapply(BootSigPath, file.exists)]
toRunBoot[BootSigExists == FALSE]
toRunBoot <- toRunBoot[BootSigExists == TRUE]



toRunBoot <- toRunBoot[,BootRes := sapply(BootSigPath, function(x) return(list(readRDS(x))))]


toRunBoot[,pvalue := sapply(BootRes, function(x) return(x[[1]]))]
toRunBoot[,estimate := sapply(BootRes, function(x) return(x[[2]]))]

toRunBoot[, lower := sapply(BootRes, function(x) return(x[["ci"]]$percent[4]))]
toRunBoot[, upper := sapply(BootRes, function(x) return(x[["ci"]]$percent[5]))]


toRunBoot[order(pvalue)]
toRunBootOut <- toRunBoot

toRunBootOut <- toRunBootOut[,BootRes:=NULL]
toRunBootOut <- toRunBootOut[,BootSigExists:=NULL]
toRunBootOut <- toRunBootOut[,BootSigPath:=NULL]
toRunBootOut <- toRunBootOut[,hetFilePath:=NULL]

toRunBoot[, FWER_genes := pmin(pvalue * corP,1)]


fwrite(toRunBootOut,"cnvRes/cnv/biomarker_res/allRes.csv")
fwrite(toRunBootOut, file="cnvRes/cnv/biomarker_res/cnv_gene_drug_tissue_res_pharmacodb.csv")
