


library(PharmacoGx)
library(data.table)
library(doParallel)

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- file.path(scratch, "Tissue_Biomarker", "rna")


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

## 86 core hours

makeToRunTable <- function(pset.list) {
	all.drugs <- .unionList(lapply(pset.list, drugNames))

	drug.table <- sapply(seq_along(pset.list), function(i) {
		return(all.drugs %in% drugNames(pset.list[[i]]))
	})

	rownames(drug.table) <- all.drugs
	colnames(drug.table) <- lapply(pset.list, name)

	myord <- order(rowSums(drug.table), decreasing = TRUE)

	drug.table <- drug.table[myord, ]
	# write.csv(drug.table, file="drugIntersectTable.csv")

	all.tissues <- .unionList(lapply(pset.list, function(x) {
		return(cellInfo(x)$tissueid)
	}))
	tissue.table <- sapply(pset.list, function(pset) {
		return(sapply(all.tissues, function(x) sum(x == cellInfo(pset)$tissueid, na.rm = TRUE)))
	})

	rownames(tissue.table) <- all.tissues
	colnames(tissue.table) <- lapply(pset.list, name)

	## using 20 cell lines as cutoff
	tissue.table <- tissue.table >= 20

	myord <- order(rowSums(tissue.table), decreasing = TRUE)
 	tissue.table <- tissue.table[myord, ]

	## Lets extract the gene expresssion values now.
	pset.list.genexp <- lapply(pset.list, function(pset) summarizeMolecularProfiles(pset, mDataNames(pset)))


	pset.list.genexp <- lapply(pset.list.genexp, function(SE) {
		rownames(SE) <- gsub(rep = "", x = rownames(SE), pat = "\\.[0-9]+$")
		return(SE) # nolint
	})
	names(pset.list.genexp) <- names(pset.list)

	pset.list.genexp <- lapply(names(pset.list.genexp), function(x) {
        ## microarray and rnaseq annotations have different column names
        gene_type_col <- ifelse("GeneBioType" %in% colnames(rowData(pset.list.genexp[[x]])), "GeneBioType", "gene_type")
        ## limiting feature space for power
        ft <- rownames(rowData(pset.list.genexp[[x]]))[rowData(pset.list.genexp[[x]])[[gene_type_col]] %in% "protein_coding"]
        return(pset.list.genexp[[x]][ft, ])
	})
    names(pset.list.genexp) <- names(pset.list)

    all.genes <- .unionList(lapply(pset.list.genexp, rownames))
    genes.table <- sapply(pset.list.genexp, function(SE) {
        return(all.genes%in% rownames(SE))
    })
    rownames(genes.table) <- all.genes
	colnames(genes.table) <- names(pset.list.genexp)

    tissues_in_3 <- rownames(tissue.table)[apply(tissue.table, 1, function(x) sum(x) >= 3)]
	drugs_in_3 <- rownames(drug.table)[apply(drug.table, 1, function(x) sum(x) >= 3)]
    genes_in_3 <- rownames(genes.table)[apply(genes.table, 1, function(x) sum(x) >= 3)]


    pset.list.genexp.m <- lapply(names(pset.list.genexp), function(x) {
        ## microarray and rnaseq annotations have different column names
        gene_type_col <- ifelse("GeneBioType" %in% colnames(rowData(pset.list.genexp[[x]])), "GeneBioType", "gene_type")
        ## limiting feature space for power
        ft <- rownames(rowData(pset.list.genexp[[x]]))[rowData(pset.list.genexp[[x]])[[gene_type_col]] %in% "protein_coding"]


        return(cbind("PSet" = x, reshape2::melt(as.is = TRUE, SummarizedExperiment::assay(pset.list.genexp[[x]])[ft, ])))
    })


    pset.genexp.dt <- rbindlist(pset.list.genexp.m)


    colnames(pset.genexp.dt)[2:3] <- c("geneid", "cellid")



    pset.genexp.dt <- pset.genexp.dt[!is.na(value), ]


    sens.num.dt <- rbindlist(lapply(pset.list, function(x) {
        cbind(
            "PSet" = name(x),
            reshape2::melt(summarizeSensitivityProfiles(x, "aac_recomputed"))
        )
    }))

    colnames(sens.num.dt) <- c("PSet", "Drug", "cellid", "value")
    sens.num.dt <- sens.num.dt[Drug %in% drugs_in_3]
    sens.num.dt.filt <- sens.num.dt[!is.na(value)]

    tissueCell.dt <- rbindlist(lapply(pset.list, function(x) {
        return(data.frame(cellNames(x), cellInfo(x)[, "tissueid"]))
    }))

    tissueCell.dt <- unique(tissueCell.dt)

    colnames(tissueCell.dt) <- c("cellid", "tissueid")

    sens.num.dt.filt <- merge(sens.num.dt.filt, tissueCell.dt, by = "cellid")   
        
    ## First, we filter to at least 20 cell lines in the tissue
    sens.num.dt.filt <- merge(sens.num.dt.filt[, .N, .(tissueid, Drug, PSet)][N >= 20], sens.num.dt.filt, by = c("tissueid", "Drug", "PSet"))
    sens.num.dt.filt[, N := NULL]

    ## now, we filter to at least 4 of cell lines with > 5% response in the tissue. 10% would be what DSS uses,
    ## but we are being a bit lenient here, since its very unlikely that artefacts from a single dataset
    ## would survive meta-analysis
    sens.num.dt.filt <- merge(sens.num.dt.filt[, sum(value > 5), .(tissueid, Drug, PSet)][V1 >= 4], sens.num.dt.filt, by = c("tissueid", "Drug", "PSet"))
    sens.num.dt.filt[, V1 := NULL]

    # Finally, we can drop the value column.
    sens.num.dt.filt[, value := NULL]

    sens.num.dt.filt[,Drug := as.character(Drug)]
    setkey(sens.num.dt.filt, "Drug","tissueid", "PSet")

    pset.genexp.dt <- merge(pset.genexp.dt, tissueCell.dt, by='cellid')

    setkey(pset.genexp.dt, geneid, tissueid, PSet)

    
    toCheckTable <- data.table(expand.grid(drugs_in_3, tissues_in_3, genes_in_3, stringsAsFactors = FALSE))

    names(toCheckTable) <- c("Drug", "Tissue", "Gene")
    

    genes.table.dt <- as.data.table(genes.table, keep.rownames = "geneid")
    gene.pset.map <- melt(genes.table.dt,
        id.vars = "geneid",
        variable.name = "PSet",
        value.name = "present",
        variable.factor = FALSE # Keep PSet as character
    )[present == TRUE, .(geneid, PSet)] # Keep only TRUE entries
    setkey(gene.pset.map, geneid) # Key for fast lookup by gene
    rm(genes.table.dt) # Clean up intermediate table

    check.row <- function(cur.drug, cur.tissue, cur.gene) {
        
        # cur.tissue.cellines <- tissueCell.dt[tissueid == cur.tissue, cellid]

        cur.sens.tbl <- sens.num.dt.filt[.(cur.drug, cur.tissue)]

        # first, we figure out which datasets this drug was tested in
        cur.drug.dataset <- cur.sens.tbl[, unique(PSet)]
        if(length(cur.drug.dataset)<3) {
            return(character(0))
        }
    
        # now, we figure out which datasets this gene was tested in
        cur.gene.dataset <- gene.pset.map[.(cur.gene), PSet, nomatch = 0L]
        if(length(cur.gene.dataset)<3) {
            return(character(0))
        }
        common.datasets <- intersect(cur.drug.dataset, cur.gene.dataset)

        if(length(common.datasets)<3) {
            return(character(0))
        }
        
        cur.genexp.dt <- pset.genexp.dt[.(cur.gene, cur.tissue)]


        common.measurements.in.dataset <- sapply(common.datasets, function(ds) {
            # dataset <- pset.list.genexp[[ds]]
            # cur.gene.values <- assay(dataset, 1)[cur.gene, intersect(cur.tissue.cellines, colnames(dataset))]
            # gene.cellines <- names(cur.gene.values)[!is.na(cur.gene.values)]
            gene.cellines <- cur.genexp.dt[PSet == ds, cellid]
            if (length(gene.cellines) < 20) {
                return(FALSE)
            }
            sens.celllines <- cur.sens.tbl[PSet == ds, cellid]
            return(length(intersect(sens.celllines, gene.cellines)) >= 20)
        })

        passing.datasets <- common.datasets[common.measurements.in.dataset]

        if(length(passing.datasets)<3) {
            return(character(0))
        }

        return(passing.datasets)
    }

    toCheckTableList <- split(toCheckTable,toCheckTable$Gene)

    registerDoParallel(6)
    out <- foreach(toCheckTable = toCheckTableList) %dopar% {
        toCheckTable[, check.row(Drug, Tissue, Gene), by = .(Drug, Tissue, Gene)]
    }

    out.table <- rbindlist(out)

    # ## this can be further optimized by execution order I think, for example, the large subsets of pset.genexp.dt by 


    # Rprof()
    # my.x <- toCheckTable[1:5000, check.row(Drug, Tissue, Gene), by = .(Drug, Tissue, Gene)]
    # Rprof(NULL)
    # summaryRprof()$by.self

    # system.time(toCheckTable[1:5000, check.row(Drug, Tissue, Gene), by = .(Drug, Tissue, Gene)])

    # p <- profvis({
    # my.x <- toCheckTable[1:10000, check.row3(Drug, Tissue, Gene), by = .(Drug, Tissue, Gene)]
    # })
    # p
    return(out.table)
}

# write.csv(tissue.table, file="tissueIntersectTable.csv")


all.dt <- makeToRunTable(pset.list)


all.dt <- all.dt[order(V1, Tissue, Drug, Gene)]

if(!dir.exists(file.path(project, "runlist_files/"))) dir.create(file.path(project, "runlist_files/"))

write.table(all.dt, file.path(project, "runlist_files/geneExpressionMasterToRunList.txt"), quote = FALSE, row.names = FALSE, sep = ",", col.names = FALSE)



