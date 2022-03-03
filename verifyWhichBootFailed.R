library(data.table)

method <- "perm"

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))

toRunBoot <- fread("~/toRunMetaByGeneBoot.txt", header=FALSE)

colnames(toRunBoot) <- c("Drug", 'Tissue', "Gene")

toRunBoot[,hetFilePath := file.path(scratch, "perm_meta_hetPermTest", 
		paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_hetPerm_10000_out.rds"))]

toRunBoot[,hetTestExists := file.exists(hetFilePath)]

toRunBoot[hetTestExists == TRUE]


# bootOut_Afatinib_Breast_ENSG00000092850_10000000_out.txt

juliaDir <- "/cluster/projects/bhklab/tmp/psmirnov_scratch/juliaBoot"

toRunBoot[,BootFilePath := file.path(juliaDir, paste0("bootOut_", make.names(Drug), "_", make.names(Tissue), "_", Gene, "_10000000_out.txt"))]




toRunBoot[,BootTestExists := sapply(BootFilePath, file.exists)]
toRunBoot[,I := .I-1]
toRerun = toRunBoot[BootTestExists != TRUE, c(1,2,3,8)]
toRerun <- data.frame(toRerun)
colnames(toRerun) <- NULL

write.csv(toRerun, "~/toRerunBoot.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

## SOME OF THESE ARE failing due to almost 0 expression in all samples in a dataset (so model fit fails)


bootSigDir <- file.path(scratch, "perm_meta_boot_sig")

toRunBoot[,BootSigPath := file.path(bootSigDir, paste0("bootSig_", make.names(Drug), "_", make.names(Tissue), "_", Gene, "_out.rds"))]

toRunBoot[,BootSigExists := sapply(BootSigPath, file.exists)]
toRunBoot[BootSigExists == FALSE]

write.csv(toRerun, "~/toRerunBootSig.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

