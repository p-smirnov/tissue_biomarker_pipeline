scratch <- Sys.getenv("SCRATCH")
res.dir <- "pearson_perm_res"

toRunEx <- read.csv("tissueDrugToRunExtended.txt", header=FALSE)



toRunEx <- do.call(function(...) paste(..., sep="_"), toRunEx)

toRunEx <- paste("signature", toRunEx, sep="_")
toRunEx <- make.names(toRunEx)


toRunEx <- paste(toRunEx, "rds", sep=".")

myfls <- list.files(file.path(scratch, res.dir))

toRunEx[which(!toRunEx %in% myfls)]

toRunAgain <- read.csv("tissueDrugToRunExtended.txt", header=FALSE)

toRunAgain <- toRunAgain[which(!toRunEx %in% myfls),]

write.table(toRunAgain, "tissueDrugToRunAgain.txt", quote=FALSE, row.names=FALSE, sep=",", col.names=FALSE)

toRunMeta <- read.csv("tissueDrugToRunExtended.txt", header=FALSE)
library(data.table)
toRunMeta <- data.table(toRunMeta)
colnames(toRunMeta) <- c("PSet", "Drug", "Tissue")
toRunMeta <- toRunMeta[,.(Drug, Tissue)]

toRunMeta <- unique(toRunMeta)


write.table(toRunMeta, "tissueDrugToRunMeta.txt", quote=FALSE, row.names=FALSE, sep=",", col.names=FALSE)
