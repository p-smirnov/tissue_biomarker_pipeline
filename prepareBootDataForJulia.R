library(PharmacoGx)
library(meta)
library(SummarizedExperiment)
library(coop)
# library(foreach)
# library(doParallel)
# library(doRNG)
library(iterators)
library(lme4)

library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
RhpcBLASctl::blas_set_num_threads(1)
# registerDoParallel(40)



args <- commandArgs(trailingOnly = TRUE)

method <- "perm"
hetTestCutoff <- 0.1
## ACSS2 is interesting



# psetName <- args[1]

# drug <- "Lapatinib"


home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myOutDir <- file.path(scratch, paste0(method, "_meta_bootstrap"))



loadPSet <- function(psetName, tissue){


	switch(psetName, 
		   CCLE = {
		   		pset <- readRDS(file.path(myDataDir,"CCLE.rds"))
		   }, CCLE.CTRPv2 = {
		   		pset <- readRDS(file.path(myDataDir,"CCLE.CTRPv2.rds"))
		   }, GDSC_v1 = {
		   		pset <- readRDS(file.path(myDataDir,"GDSC1.rds"))
		   }, GDSC_v2 = {
		   		pset <- readRDS(file.path(myDataDir,"GDSC2.rds"))
		   }, gCSI = {
		   		pset <- readRDS(file.path(myDataDir,"gCSI.rds"))
		   }, GRAY = {
		   		pset <- readRDS(file.path(myDataDir,"GRAY.rds"))
		   }, UHNBreast = {
		   		pset <- readRDS(file.path(myDataDir,"UHNBreast.rds"))
		   }, {stop("Please Provide a valid pset")})

	mData <- mDataNames(pset)


	gene_type_col <- ifelse("GeneBioType" %in% colnames(featureInfo(pset, mData)), "GeneBioType", "gene_type")

	ft <- rownames(featureInfo(pset, mData))[which(featureInfo(pset, mData)[[gene_type_col]] == "protein_coding")]

	if(is.na(tissue)){
		chosen.cells <- cellNames(pset)
	} else {
		chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
	}

	return(list(pset = pset, mData = mData, chosen.cells = chosen.cells))
}



getTrange <- function(df){
	return(qt(0.975, df) - qt(0.025, df))
}


toRunExtended <- read.csv("~/tissueDrugToRunExtended.txt", header=FALSE)
pSets <- toRunExtended[toRunExtended[,2] == drug & toRunExtended[,3] == tissue, 1]

geneInfo <- read.csv("~/geneInfo.csv")


if(!file.exists(myOutDir)) dir.create(myOutDir)

# fls <- list.files(path=), pattern = paste0("*_", make.names(drug), "_", make.names(tissue), "*"), full.names = TRUE)
# fls <- paste0(file.path(mySigDir,paste("signature", make.names(pSets), make.names(drug), make.names(tissue), sep="_")), ".rds")
# sig.list <- sapply(fls, readRDS)

# genes_to_check <- unique(unlist(sapply(sig.list, function(x) {
# 	mygenes <- names(which(1==x[,1,"significant"]))
# 	mygenes <- gsub(mygenes, pat="\\.[0-9]+$", rep="")
# 	return(mygenes)
# 	})))


# psetsRequired <- sapply(sig.list, function(x) return(x@PSetName))

pset.list <- lapply(pSets, loadPSet, tissue=tissue)
names(pset.list) <- pSets


getBootSample <- function(model.data, scale=TRUE){

	sampled.datasets <- sample(unique(model.data$dataset), replace=TRUE)
	
	names(sampled.datasets) <- make.unique(sampled.datasets)

	sampled.data <- lapply(names(sampled.datasets), function(dt){
		sm.dt <- model.data[sample(which(model.data$dataset == sampled.datasets[dt]), replace=TRUE), ]
		sm.dt$dataset <- dt
        sm.dt[,1] <- scale(sm.dt[,1])
        sm.dt[,2] <- scale(sm.dt[,2])
		return(sm.dt)
		})
	sampled.data <- do.call(rbind, sampled.data)

	return(sampled.data)
}

bootMeanEffects <- function(model.data,R){

    N <- nrow(model.data)

    t <- unlist(mclapply(seq_len(R), function(i, model.data) {
    	library(lme4)
        test <- getBootSample(model.data)
        fixedEff <- tryCatch({
            test.mod <- lmer(y ~ (x + 0| dataset) + x + 0, test, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            fixef(test.mod)

            }, error = function(e){ NA_real_})
        fixedEff
    }, model.data = model.data, mc.cores=nthread))


    ## Calculating t0
    model.data.scaled <- do.call(rbind, lapply(model.data$dataset, function(dt){
        sm.dt <- model.data[model.data$dataset == dt,]
        sm.dt[,1] <- scale(sm.dt[,1])
        sm.dt[,2] <- scale(sm.dt[,2])
        return(sm.dt)
        }))


    t0 <- fixef(lmer(y ~ (x + 0| dataset) + x + 0, model.data.scaled))



	dim(t) <- c(R,1)
	res <- list(
	    t0 = t0,
	    t = t,
	    R = R,
	    data = model.data,
	    seed = seed,
	    sim = "ordinary",
	    stype = "i",
	    call = match.call(),
	    strata = as.numeric(as.factor(model.data$dataset)),
	    weights = rep(1/N, N)
	)
    attr(res, "class") <- "boot"
    attr(res, "boot_type") <- "boot"
    return(res)
}

# gene <- "ENSG00000141736"


permFixedEffect <- function(model.data, R){

	# t <- foreach(i =icount(R), .combine=c, .inorder=FALSE) %dorng% {
	# 	library(RhpcBLASctl)
	# 	RhpcBLASctl::blas_set_num_threads(1)
 #        x <- model.data[,"x"]
 #        y <- sample(model.data[,"y"])
 #        # coop::pcor(x,y)
 #        cor(x,y)
 #    }
	t <- unlist(mclapply(seq_len(R), function(i, model.data) {
        x <- model.data[,"x"]
        y <- sample(model.data[,"y"])
        # coop::pcor(x,y)
        cor(x,y)
	}, model.data = model.data, mc.cores=nthread))
    x <- model.data[,"x"]
	y <- model.data[,"y"]
    t0 <- coop::pcor(x,y)
    # t0 <- cor(x,y)
    return(list(t0 = t0, t = t, R = R, data = model.data))
}

x.list <- lapply(pset.list, function(pset.pars){

	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	mData <- pset.pars$mData

	mol.prof <- assay(summarizeMolecularProfiles(pset, mData, cell.lines = chosen.cells))

	myx <- grep(gene, rownames(mol.prof))

	if(!length(myx)) return(rep(NA_real_, times=length(pset.pars$chosen.cells)))
	if(length(myx) > 1) stop(paste0("multiple genes matched:", gene, " id. please investigate"))

	mol.prof <- mol.prof[myx,]

	return(mol.prof)


	})

y.list <- lapply(pset.list, function(pset.pars){

	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	
	drug.res <- summarizeSensitivityProfiles(pset, "aac_recomputed", cell.lines = chosen.cells)[drug,]


	return(drug.res)


	})

model.data <- data.frame(x = unlist(x.list),
						 y = unlist(y.list),
						 dataset = unlist(lapply(pset.list, function(pset.pars){
						 	return(rep(name(pset.pars$pset), times=length(pset.pars$chosen.cells)))
						 	})))


model.data <- model.data[complete.cases(model.data),]



gene.fls <- paste0("~/featureLists/", make.names(tissue), ".", make.names(pSets), ".csv")
total_gene_list <- .unionList(lapply(gene.fls, readLines))
corrected_alpha <- corrected_alpha/length(total_gene_list)

R <- ceiling(1/corrected_alpha)
rm(pset.list, y.list, x.list)
gc()



