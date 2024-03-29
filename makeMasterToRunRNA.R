


library(PharmacoGx)


home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "TBPInputs", "rna")

CCLE <- readRDS(file.path(myDataDir,"CCLE.rds"))
CTRPv2 <- readRDS(file.path(myDataDir,"CCLE.CTRPv2.rds"))
GDSC1 <- readRDS(file.path(myDataDir,"GDSC1.rds"))
GDSC2 <- readRDS(file.path(myDataDir,"GDSC2.rds"))
gCSI <- readRDS(file.path(myDataDir,"gCSI.rds"))
GRAY <- readRDS(file.path(myDataDir,"GRAY.rds"))
UHNBreast <- readRDS(file.path(myDataDir,"UHNBreast.rds"))

pset.list <- list(CCLE, CTRPv2, GDSC1, GDSC2, GRAY, UHNBreast, gCSI)

names(pset.list) <- sapply(pset.list, name)

all.drugs <- .unionList(lapply(pset.list, drugNames))


drug.table <- sapply(seq_along(pset.list), function(i) return(all.drugs %in% drugNames(pset.list[[i]])))

rownames(drug.table) <- all.drugs
colnames(drug.table) <- lapply(pset.list, name)

myord <- order(rowSums(drug.table), decreasing=TRUE)

drug.table <- drug.table[myord,]
write.csv(drug.table, file="drugIntersectTable.csv")

all.tissues <- .unionList(lapply(pset.list, function(x) return(cellInfo(x)$tissueid)))
tissue.table <- sapply(pset.list, function(pset) return(sapply(all.tissues, function(x) sum(x == cellInfo(pset)$tissueid, na.rm=TRUE))))

rownames(tissue.table) <- all.tissues
colnames(tissue.table) <- lapply(pset.list, name)

## using 20 cell lines as cutoff
tissue.table <- tissue.table >= 20

myord <- order(rowSums(tissue.table), decreasing=TRUE)
tissue.table <- tissue.table[myord,]
write.csv(tissue.table, file="tissueIntersectTable.csv")


library(reshape2)
library(data.table)
library(SummarizedExperiment)
## Final file output here


## TODO: run the full damn merge and see if the numbers match

### Approximate design of this code: filter to drug intersection in at least 3 datasets
### Filter to tissues with >20 cell lines in at least 3 datasets
### 


tissues_in_3 <- rownames(tissue.table)[apply(tissue.table, 1, function(x) sum(x)>=3)]
drugs_in_3 <- rownames(drug.table)[apply(drug.table, 1, function(x) sum(x)>=3)]

cells_per_dataset <- sapply(pset.list, function(pset) return(cellNames(pset)[cellInfo(pset)$tissueid %in% tissues_in_3]))

sens.num.dt <- rbindlist(lapply(pset.list, function(x) cbind("PSet" = name(x), 
	reshape2::melt(summarizeSensitivityProfiles(x, "aac_recomputed")))))


colnames(sens.num.dt) <- c("PSet", "Drug", "cellid", "value")

sens.num.dt <- sens.num.dt[cellid %in% .unionList(cells_per_dataset)]
sens.num.dt <- sens.num.dt[Drug %in% drugs_in_3]

sens.num.dt <- sens.num.dt[!is.na(value)]

## Here we filtered to cell lines in tissues that are in 3+ datasets, and drugs tested on these cell. 
## We drop the value column, and keep only the cells-drug combinations that occur at least 3 times

sens.num.dt[,value := NULL]

sens.num.dt.filt <- sens.num.dt
## TODO:: This is wrong, I don't actually want the cellid to be in more than three datasets, just more than 20 cells in at least 3 datasets
# sens.num.dt.filt <- merge(sens.num.dt[,.N, .(cellid, Drug)][N>=3,.(cellid, Drug)], sens.num.dt)

## Now lets merge this with tissue information for each cell line

tissueCell.dt <- rbindlist(lapply(pset.list, function(x) return(data.frame(cellNames(x), cellInfo(x)[,"tissueid"]))))

tissueCell.dt <- unique(tissueCell.dt)

colnames(tissueCell.dt) <- c("cellid", "tissueid")

sens.num.dt.filt <- merge(sens.num.dt.filt, tissueCell.dt, by="cellid")

sens.num.dt.filt <- merge(sens.num.dt.filt[,.N,.(tissueid, Drug, PSet)][N>=20], sens.num.dt.filt, by=c("tissueid", "Drug", "PSet"))
sens.num.dt.filt[,N:=NULL]
## Lets extract the gene expresssion values now.
pset.list.genexp <- lapply(pset.list, function(pset) summarizeMolecularProfiles(pset, mDataNames(pset)))


pset.list.genexp <- lapply(pset.list.genexp, function(SE) {
	rownames(SE) <- gsub(rep="", x=rownames(SE), pat="\\.[0-9]+$")
	return(SE)
})
names(pset.list.genexp) <- names(pset.list)




pset.list.genexp.m <- lapply(names(pset.list.genexp), function(x) return(cbind("PSet" = x,reshape2::melt(as.is =TRUE, SummarizedExperiment::assay(pset.list.genexp[[x]])))))
pset.genexp.dt <- rbindlist(pset.list.genexp.m)


colnames(pset.genexp.dt)[2:3] <- c("geneid", "cellid")



pset.genexp.dt <- pset.genexp.dt[!is.na(value),]

setkey(pset.genexp.dt, geneid)

## I will now filter out geneids not occuring in at least 3 datasets


pset.genexp.dt <- merge(pset.genexp.dt,pset.genexp.dt[,length(unique(PSet)), geneid][V1>=3, .(geneid)])


## Now filter to only cell lines in my tissues of interest

pset.genexp.dt <- pset.genexp.dt[cellid %in% unique(sens.num.dt.filt$cellid)]

cellTissue.dt <- unique(sens.num.dt.filt[,.(cellid, tissueid)])

pset.genexp.dt <- merge(pset.genexp.dt, cellTissue.dt, by="cellid")


## Now check that per-dataset, we have 20 cell lines per tissue per gene



pset.genexp.dt <- merge(pset.genexp.dt, pset.genexp.dt[,length(unique(cellid)), .(PSet, tissueid, geneid)][V1 >= 20], by=c("PSet", "tissueid", "geneid"))



pset.genexp.dt[,V1 := NULL]

setkey(pset.genexp.dt, cellid)

pset.genexp.dt[,value := NULL]

## Need to do this by tissue, as I run out of memory


pset.genexp.dt.split <- split(pset.genexp.dt, by="tissueid")

sens.num.dt.filt.split <- split(sens.num.dt.filt, by="tissueid")

all.dt.list <- list()

for(tissue in names(pset.genexp.dt.split)){


	pset.genexp.dt.split[[tissue]]
	all.dt.tissue <-merge(pset.genexp.dt.split[[tissue]], sens.num.dt.filt.split[[tissue]], by=c("cellid", "PSet", "tissueid"),allow.cartesian=TRUE)

	all.dt.tissue <- merge(all.dt.tissue, all.dt.tissue[,.N, .(Drug, geneid, tissueid, PSet)][N>=20,], by=c("Drug", "tissueid", "geneid", "PSet"))
	### now that I am sure there are 20 samples per pset for each drug, gene, tissue triplet, I can drop the cellid column

	all.dt.tissue[,cellid:=NULL]

	all.dt.tissue[,N:=NULL]

	all.dt.tissue <- unique(all.dt.tissue)

	all.dt.tissue <- merge(all.dt.tissue, all.dt.tissue[,.N, .(geneid, tissueid, Drug)][N >= 3], by=c("geneid", "tissueid", "Drug"))

	all.dt.list[[tissue]] <- all.dt.tissue
	gc()
}


all.dt <- rbindlist(all.dt.list)


all.dt[,N:=NULL]
all.dt <- all.dt[order(PSet, tissueid, Drug, geneid)]

write.table(all.dt, file.path(scratch, "geneExpressionMasterToRunList.txt"), quote=FALSE, row.names=FALSE, sep=",", col.names=FALSE)



