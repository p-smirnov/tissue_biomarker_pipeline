library(data.table)

method <- "perm"



home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "smallPSets")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_perm_sig"))
myPvalDir <- file.path("/zfs-data/tbp_prism_res", paste0(method, "_meta_perm_sig"))
badchars <- "[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[ ]|[(]|[)]"

make.names.2 <- function(x) {
    return(gsub(pat = badchars, rep = ".", x))
}

torunall <- fread("/zfs-data/tbp_prism_res/runlist_files/geneExpressionMasterToRunList.txt",header=F)

corP <- length(unique(torunall[,V1]))

## First we load in the fixed effect permutation test results:

toRunAfterHet <- fread("/zfs-data/tbp_prism_res/runlist_files/metaHetTestRes.txt")

toRunPerm <- toRunAfterHet[hetTestRes==0]

# colnames(toRunPerm) <- c("Drug", 'Tissue', "Gene")

# toRunPerm[,hetFilePath := file.path(scratch, "perm_meta_hetPermTest", 
# 		paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_hetPerm_10000_out.rds"))]

# toRunPerm[,hetTestExists := file.exists(hetFilePath)]

# toRunPerm[hetTestExists == TRUE]


# toRunPerm[,PermFilePattern := paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_perm*")]

# toRunPerm[,PermFilePath := sapply(PermFilePattern, function(x) list.files(path=file.path(scratch, "perm_meta_bootstrap"), pattern=x, full.names=TRUE))]

# toRunPerm[,PermTestExists := sapply(PermFilePath, file.exists)]
# toRunPerm[!sapply(PermTestExists, any)]

toRunPerm[, PermResFile := file.path(myPvalDir, paste0("permSig_", Gene, "_", make.names.2(Drug), "_", make.names.2(Tissue), "_out.rds"))]
toRunPerm[,all(sapply(PermResFile, file.exists))]

system.time(toRunPerm[,PermRes := sapply(PermResFile, function(x) list(readRDS(x)))])

toRunPerm[,pvalue := sapply(PermRes, function(x) return(x[[1]]))]
toRunPerm[,estimate := sapply(PermRes, function(x) return(x[[2]]$t0))]

toRunPerm[, lower := sapply(PermRes, function(x) return(x[["ci"]]$perc[4]))]
toRunPerm[, upper := sapply(PermRes, function(x) return(x[["ci"]]$perc[5]))]


toRunPerm[order(pvalue)]

toRunPermOut <- toRunPerm

toRunPermOut <- toRunPermOut[,PermRes:=NULL]
toRunPermOut <- toRunPermOut[,PermResFile:=NULL]

toRunPermOut <- toRunPermOut[, hetTestRes := NULL]
toRunPermOut <- toRunPermOut[, hetFilePath := NULL]

write.csv(toRunPermOut, file = "/zfs-data/tbp_prism_res/biomarker_res/BiomarkerPermRes.csv")


### Now lets do the bootstrap results



# home <- Sys.getenv("HOME")
# scratch <- Sys.getenv("SCRATCH")

# myDataDir <- file.path(home, "Data", "smallPSets")
# mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
# bootSigDir <- file.path(scratch, "perm_meta_boot_sig")
bootSigDir <- file.path("/zfs-data/tbp_prism_res", paste0(method, "_meta_boot_sig"))



# toRunBoot <- fread("toRunFiles/toRunMetaByGeneBoot.txt", header=FALSE)

# toRunAfterHet <- fread("/zfs-data/tbp_prism_res/runlist_files/metaHetTestRes.txt")

toRunBoot <- toRunAfterHet[hetTestRes == 1]

# colnames(toRunBoot) <- c("Drug", 'Tissue', "Gene")




toRunBoot[, BootSigPath := file.path(bootSigDir, paste0("bootSig_", Gene, "_", make.names.2(Drug), "_", make.names.2(Tissue), "_out.rds"))]

toRunBoot[,BootSigExists := sapply(BootSigPath, file.exists)]
toRunBoot[BootSigExists == FALSE]
toRunBoot <- toRunBoot[BootSigExists == TRUE]
## some jobs seem to have failed, run them

toRunBoot <- toRunBoot[,BootRes := sapply(BootSigPath, function(x) return(list(readRDS(x))))]


toRunBoot[,pvalue := sapply(BootRes, function(x) return(x[[1]]))]
toRunBoot[,estimate := sapply(BootRes, function(x) return(x[[2]]))]

toRunBoot[, lower := sapply(BootRes, function(x) return(x[["ci"]]$perc[4]))]
toRunBoot[, upper := sapply(BootRes, function(x) return(x[["ci"]]$perc[5]))]


toRunBoot[order(pvalue)]
toRunBootOut <- toRunBoot

toRunBootOut <- toRunBootOut[,BootRes:=NULL]
toRunBootOut <- toRunBootOut[,BootSigExists:=NULL]
toRunBootOut <- toRunBootOut[,BootSigPath:=NULL]

toRunBootOut <- toRunBootOut[, hetTestRes := NULL]
toRunBootOut <- toRunBootOut[, hetFilePath := NULL]

write.csv(toRunBootOut, file = "/zfs-data/tbp_prism_res/biomarker_res/BiomarkerBootRes.csv")

all.res <- rbind(toRunBootOut,toRunPermOut)

all.res[, FWER_genes := pmin(pvalue * corP,1)]


gene_info <- fread("geneInfo.csv")

gene_info[, V1 := gsub(V1, pat = "\\.[0-9]+$", rep = "")]

gene_symbol <- gene_info[, .(V1, gene_name)]

colnames(gene_symbol)[1] <- "Gene"

all.res <- merge(all.res, gene_symbol, by = "Gene", all.x = TRUE)

fwrite(all.res, file = "rnaResultsPrism/biomarker_res/mrna_gene_drug_tissue_res_pharmacodb.csv")




