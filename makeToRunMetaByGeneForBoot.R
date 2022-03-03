library(PharmacoGx)
library(data.table)

method <- "perm"
args <- commandArgs(trailingOnly = TRUE)
outDir <- args[[2]]
inDir <- args[[1]]

home <- Sys.getenv("HOME")
project <- Sys.getenv("PROJECT")


badchars <- "[ ]|[/]|[:]|[-]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))



myToRunFileName <- file.path(inDir, "toRunMetaByGene.txt")

hetTestDir <- file.path(project, paste0("hetTest"))


## dealing with possibly running in the cloud
containername <- Sys.getenv("containername", unset=NA_character_)

if(!is.na(containername)){
	hetTestDir <- file.path(containername, hetTestDir)
	inDir <- file.path(containername, inDir)
	outDir <- file.path(containername, outDir)

	myToRunFileName <- file.path(containername, myToRunFileName)
}


toRun <- fread(myToRunFileName, header=FALSE)

colnames(toRun) <- c("Gene", "Tissue", "Drug", "PSet", "Significant")
toRunHet <- unique(toRun[Significant==1,.(Gene, Tissue, Drug)])

 # TODO: make the R for het testing a parameter that can be passed in 
toRunHet[,hetFilePath := file.path(hetTestDir, 
		make.names.2(paste0(Drug, "_", Tissue, "_", Gene, "_hetPerm_10000_out.rds")))]

toRunHet[,hetTestExists := file.exists(hetFilePath)]

toRunMetaFinal <- toRunHet#[(hetTestExists)] #should exist for each gene!


toRunMetaFinal[,hetTestRes := sapply(hetFilePath, function(x) return(as.numeric(readRDS(x)$pvalue < 0.1)))]

toRunMetaFinal[,hetTestExists := NULL]

fwrite(toRunMetaFinal, file=file.path(outDir, "metaHetTestRes.txt"))

# toRunMetaPerm <- toRunMetaFinal[hetTestRes == FALSE]
# toRunMetaPerm <- toRunMetaPerm[,c("Drug", "Tissue", "Gene"), with=FALSE]


# fwrite(toRunMetaPerm, file="~/toRunMetaByGenePerm.txt", sep=",", col.names=FALSE, quote=FALSE)



# toRunMetaBoot <- toRunMetaFinal[hetTestRes == TRUE]
# toRunMetaBoot <- toRunMetaBoot[,c("Drug", "Tissue", "Gene"), with=FALSE]

# fwrite(toRunMetaBoot, file="~/toRunMetaByGeneBoot.txt", sep=",", col.names=FALSE, quote=FALSE)



# toRunMetaFinal[,.N, c("Drug", "Tissue")]