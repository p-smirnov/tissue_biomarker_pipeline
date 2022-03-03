

library(PharmacoGx)


home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "TBPInputs", "cnv")


CTRPv2 <- readRDS(file.path(myDataDir,"CCLE.CTRPv2.rds"))

GDSC2 <- readRDS(file.path(myDataDir,"GDSC2.rds"))
gCSI <- readRDS(file.path(myDataDir,"gCSI.rds"))



pset.list <- list(CTRPv2, GDSC2,gCSI)

names(pset.list) <- sapply(pset.list, name)

all.drugs <- .unionList(lapply(pset.list, drugNames))


drug.table <- sapply(seq_along(pset.list), function(i) return(all.drugs %in% drugNames(pset.list[[i]])))

rownames(drug.table) <- all.drugs
colnames(drug.table) <- lapply(pset.list, name)

myord <- order(rowSums(drug.table), decreasing=TRUE)

drug.table <- drug.table[myord,]



library(reshape2)
library(data.table)
library(SummarizedExperiment)
## Final file output here


## TODO: run the full damn merge and see if the numbers match

### Approximate design of this code: filter to drug intersection in at least 3 datasets
### 


drugs_in_3 <- rownames(drug.table)[apply(drug.table, 1, function(x) sum(x)>=3)]

cells_per_dataset <- sapply(pset.list, function(pset) return(cellNames(pset)))

sens.num.dt <- rbindlist(lapply(pset.list, function(x) cbind("PSet" = name(x), 
	reshape2::melt(summarizeSensitivityProfiles(x, "aac_recomputed")))))


colnames(sens.num.dt) <- c("PSet", "Drug", "cellid", "value")

sens.num.dt <- sens.num.dt[cellid %in% .unionList(cells_per_dataset)]
sens.num.dt <- sens.num.dt[Drug %in% drugs_in_3]

sens.num.dt <- sens.num.dt[!is.na(value)]

## Here we filtered to cell lines in 3+ datasets, and drugs tested on these cell. 
## We drop the value column, and keep only the cells-drug combinations that occur at least 3 times

sens.num.dt[,value := NULL]

sens.num.dt.filt <- sens.num.dt

sens.num.dt.filt <- merge(sens.num.dt.filt[,.N,.(Drug, PSet)][N>=20], sens.num.dt.filt, by=c("Drug", "PSet"))
sens.num.dt.filt[,N:=NULL]



## Load in the CGc genes

cgc.table <- fread("~/Code/tissue_biomarker/cancerGeneCensus/cgc.genes.csv")



## Lets extract the gene expresssion values now.
pset.list.genexp <- lapply(pset.list, function(pset) {
	# mData <- mDataNames(pset)
	# gene_type_col <- ifelse("GeneBioType" %in% colnames(featureInfo(pset, mData)), "GeneBioType", "gene_type") 
	
	# ## limiting feature space for power
	# ft <- rownames(featureInfo(pset, mData))[featureInfo(pset, mData)[[gene_type_col]] == "protein_coding"]
	ft <- intersect(cgc.table$Ensembl.id, rownames(featureInfo(pset, mDataNames(pset))))
	return(summarizeMolecularProfiles(pset, mDataNames(pset), features=ft))
})


# pset.list.genexp <- lapply(pset.list.genexp, function(SE) {
# 	rownames(SE) <- gsub(rep="", x=rownames(SE), pat="\\.[0-9]+$")
# 	return(SE)
# })
names(pset.list.genexp) <- names(pset.list)




pset.list.genexp.m <- lapply(names(pset.list.genexp), function(x){ 

	return(cbind("PSet" = x,reshape2::melt(as.is =TRUE, SummarizedExperiment::assay(pset.list.genexp[[x]]))))
})
pset.genexp.dt <- rbindlist(pset.list.genexp.m)


colnames(pset.genexp.dt)[2:3] <- c("geneid", "cellid")



pset.genexp.dt <- pset.genexp.dt[!is.na(value),]

setkey(pset.genexp.dt, geneid)

## I will now filter out geneids not occuring in at least 3 datasets


pset.genexp.dt <- merge(pset.genexp.dt,pset.genexp.dt[,length(unique(PSet)), geneid][V1>=3, .(geneid)])


## Now filter to only cell lines in my tissues of interest

pset.genexp.dt <- pset.genexp.dt[cellid %in% unique(sens.num.dt.filt$cellid)]


## Now check that per-dataset, we have 20 cell lines  per gene



pset.genexp.dt <- merge(pset.genexp.dt, pset.genexp.dt[,length(unique(cellid)), .(PSet, geneid)][V1 >= 20], by=c("PSet", "geneid"))



pset.genexp.dt[,V1 := NULL]

setkey(pset.genexp.dt, cellid)

pset.genexp.dt[,value := NULL]


all.dt.list <- list()

## the cartisian join and then filter is not possible, need to loop.

for(drug in unique(sens.num.dt.filt$Drug)){


	drug.sens.num <- sens.num.dt.filt[Drug == drug]
	all.dt.drug <- merge(pset.genexp.dt, drug.sens.num, by=c("cellid", "PSet"), allow.cartesian=TRUE)

	all.dt.drug <- merge(all.dt.drug, all.dt.drug[,.N, .( geneid, PSet)][N>=20,], by=c("geneid", "PSet"))


	all.dt.drug[,cellid:=NULL]

	all.dt.drug[,N:=NULL]
	all.dt.drug <- unique(all.dt.drug)

	all.dt.drug <- merge(all.dt.drug, all.dt.drug[,.N, .(geneid, Drug)][N >= 3], by=c("geneid", "Drug"))

	all.dt.list[[drug]] <- all.dt.drug

	gc()

}



all.dt <- rbindlist(all.dt.list)


all.dt[,N:=NULL]
all.dt[,tissueid:="all"]


all.dt <- all.dt[order(PSet, Drug, geneid)]

all.dt <- all.dt[,.(geneid, tissueid, Drug, PSet)]


write.table(all.dt, file.path(scratch, "cnvCGStoRunList.txt"), quote=FALSE, row.names=FALSE, sep=",", col.names=FALSE)



