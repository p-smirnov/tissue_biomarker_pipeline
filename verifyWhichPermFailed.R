library(data.table)

method <- "perm"

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_perm_sig"))

toRunPerm <- fread("~/toRunMetaByGenePerm.txt", header=FALSE)

colnames(toRunPerm) <- c("Drug", 'Tissue', "Gene")

toRunPerm[,hetFilePath := file.path(scratch, "perm_meta_hetPermTest", 
		paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_hetPerm_10000_out.rds"))]

toRunPerm[,hetTestExists := file.exists(hetFilePath)]

toRunPerm[hetTestExists == TRUE]


toRunPerm[,PermFilePattern := paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_perm*")]

toRunPerm[,PermFilePath := sapply(PermFilePattern, function(x) list.files(path=file.path(scratch, "perm_meta_bootstrap"), pattern=x, full.names=TRUE))]

toRunPerm[,PermTestExists := sapply(PermFilePath, file.exists)]
toRunPerm[!sapply(PermTestExists, any)]

## Now we verify if the plots worked 

toRunPerm[,PermResFile := file.path(myPvalDir, paste0("permSig_", make.names(Drug), "_", make.names(Tissue), "_", Gene,"_out.rds"))]
toRunPerm[,all(sapply(PermResFile, file.exists))]



