library(data.table)




method <- "perm"

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "smallPSets")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))


toRunPerm <- fread("~/toRunMetaByGenePerm.txt", header=FALSE)

colnames(toRunPerm) <- c("Drug", 'Tissue', "Gene")

toRunPerm[,PermFilePattern := paste0(make.names(Drug), "_", make.names(Tissue), "_", Gene, "_perm*")]

toRunPerm[,PermFilePath := sapply(PermFilePattern, function(x) list.files(path=file.path(scratch, "perm_meta_bootstrap"), pattern=x, full.names=TRUE))]

toRunPerm[sapply(PermFilePath, length)>1]

p.vals <- numeric(nrow(toRunPerm))

for(i in seq_len(nrow(toRunPerm))){

	test <- readRDS(toRunPerm[i,PermFilePath])
	p.vals[i] <- mean(abs(test$t0) <= abs(test$t))
	if(!i%%100) print(i)
	rm(test)
}

# toRunPerm[,PermRes := sapply(PermFilePath, readRDS)]

