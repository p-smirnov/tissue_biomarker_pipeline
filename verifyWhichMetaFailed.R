library(data.table)

method <- "perm"

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_perm_sig"))


badchars <- "[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[ ]|[(]|[)]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))


test <- fread("/scratch/b/bhaibeka/psmirnov/Tissue_Biomarker/rna/runlist_files/metaHetTestRes.txt")

myfilepath <- "/scratch/b/bhaibeka/psmirnov/Tissue_Biomarker/rna/perm_meta_perm_sig"

##  checking fixed effect res
test[hetTestRes==0, file.exists(file.path(myfilepath, paste0("permSig_", Gene, "_",make.names.2(Drug), "_", make.names.2(Tissue), "_out.rds")))]


## checking julia data output
myfilepath <- "/scratch/b/bhaibeka/psmirnov/rna/data4juliaBoot"

which(!test[hetTestRes==1, file.exists(file.path(myfilepath, paste0("modelData_", Gene, "_",make.names.2(Drug), "_", make.names.2(Tissue), ".csv")))])
