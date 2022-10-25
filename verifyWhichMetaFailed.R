library(data.table)

method <- "perm"

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_perm_sig"))


badchars <- "[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[ ]|[(]|[)]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))




test[hetTestRes==0, file.exists(file.path(myfilepath, paste0("permSig_", Gene, "_",make.names.2(Drug), "_", make.names.2(Tissue), "_out.rds")))]

