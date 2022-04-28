library(PharmacoGx)
library(meta)
library(SummarizedExperiment)
library(coop)
# library(foreach)
# library(doParallel)
# library(doRNG)
library(iterators)
library(lme4)
library(data.table)

library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
RhpcBLASctl::blas_set_num_threads(1)
# registerDoParallel(40)
dyn.load("~/tissue_biomarker/rnaPipeline/metaPermC.so")
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

method <- "perm"
hetTestCutoff <- 0.1
## ACSS2 is interesting

drug <- "Crizotinib"
tissue <- "Lymphoid"
gene <- "ENSG00000069431"

# psetName <- args[1]

# drug <- "Lapatinib"

drug <- args[1]

# tissue <- "Breast"
tissue <- args[2]

gene <- args[3]


## 206 drugs/tissues being tested + 2 orders of magnitude for accurate computations
corrected_alpha <- as.numeric(args[4])/100 

nthread <- as.numeric(args[5]) #40 threads faster than 80 on niagara

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")

myDataDir <-  Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myOutDir <- file.path(project, paste0(method, "_meta_out"))

runDir <- args[6]

badchars <- "[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[ ]|[(]|[)]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))


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


toRunMetaRes <- fread(file=file.path(runDir, "metaHetTestRes.txt"))

toRunByGene <- toRunMetaRes[hetTestRes == FALSE]


myToRunFileName <- file.path(runDir, "toRunMetaByGene.txt")

toRunExtended <- fread(myToRunFileName, header=FALSE)

colnames(toRunExtended) <- c("Gene", "Tissue", "Drug", "PSet", "Significant")

toRunMetaRes <- merge(toRunByGene, toRunExtended, on=c("Gene", "Tissue", "Drug"))
total_gene_list <- toRunExtended[,unique(Gene)]


toRunThis <- toRunMetaRes[make.names.2(Drug) == drug & make.names.2(Tissue) == tissue & make.names.2(Gene) == gene, ]

pSets <- toRunThis[,PSet]
drug <- unique(toRunThis[,Drug])
tissue <- unique(toRunThis[,Tissue])
gene <- unique(toRunThis[,Gene])

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

# bootMeanEffects <- function(model.data,R){

#     N <- nrow(model.data)

#     t <- unlist(mclapply(seq_len(R), function(i, model.data) {
#     	library(lme4)
#         test <- getBootSample(model.data)
#         fixedEff <- tryCatch({
#             test.mod <- lmer(y ~ (x + 0| dataset) + x + 0, test, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
#             fixef(test.mod)

#             }, error = function(e){ NA_real_})
#         fixedEff
#     }, model.data = model.data, mc.cores=nthread))


#     ## Calculating t0
#     model.data.scaled <- do.call(rbind, lapply(model.data$dataset, function(dt){
#         sm.dt <- model.data[model.data$dataset == dt,]
#         sm.dt[,1] <- scale(sm.dt[,1])
#         sm.dt[,2] <- scale(sm.dt[,2])
#         return(sm.dt)
#         }))


#     t0 <- fixef(lmer(y ~ (x + 0| dataset) + x + 0, model.data.scaled))



# 	dim(t) <- c(R,1)
# 	res <- list(
# 	    t0 = t0,
# 	    t = t,
# 	    R = R,
# 	    data = model.data,
# 	    seed = seed,
# 	    sim = "ordinary",
# 	    stype = "i",
# 	    call = match.call(),
# 	    strata = as.numeric(as.factor(model.data$dataset)),
# 	    weights = rep(1/N, N)
# 	)
#     attr(res, "class") <- "boot"
#     attr(res, "boot_type") <- "boot"
#     return(res)
# }

# # gene <- "ENSG00000141736"


standardizeByDatasetInPerm <- function(iny,ds.vec){
	# y <- numeric(length(iny))
	for(myx in ds.vec){
		# model.data[myx,"x"] <- scale(model.data[myx,"x"])
		iny[myx] <- scale(iny[myx])
	}
	return(iny)
}

# # # Tested implementations below compared to this simpler one using ks test on 10k samples. 
# # # KS test found no significant differences.
# permFixedEffectSimple <- function(model.data, R){

# 	datasets <- unique(model.data[,"dataset"])
	
 
# 	t <- unlist(mclapply(seq_len(R), function(i, model.data, datasets) {
#         # model.data[,"x"] <- model.data[,"x"]
#         model.data[,"y"] <- sample(model.data[,"y"])

# 		for(ds in datasets){
# 			myx <- model.data[,"dataset"] == ds
# 			model.data[myx,"x"] <- scale(model.data[myx,"x"],scale=F)
# 			model.data[myx,"y"] <- scale(model.data[myx,"y"],scale=F)	

# 		}
# 		return(cor(model.data[,"x"], model.data[,"y"]))
# 	}, model.data = model.data, datasets = datasets, mc.cores=nthread))

#      t0 <- NA_real_
#     # t0 <- cor(x,y)
#     return(list(t0 = t0, t = t))
# }


# # # ## Current timing, 7 seconds per 5e3, N=320
# permFixedEffect <- function(model.data, R){

# 	# t <- foreach(i =icount(R), .combine=c, .inorder=FALSE) %dorng% {
# 	# 	library(RhpcBLASctl)
# 	# 	RhpcBLASctl::blas_set_num_threads(1)
#  #        x <- model.data[,"x"]
#  #        y <- sample(model.data[,"y"])
#  #        # coop::pcor(x,y)
#  #        cor(x,y)
#  #    }
# 	datasets <- unique(model.data[,"dataset"])
# 	for(ds in datasets){
# 		myx <- model.data[,"dataset"] == ds
# 		model.data[myx,"x"] <- scale(model.data[myx,"x"])
# 	}
# 	ds.vec <- lapply(datasets, function(ds) return(model.data[,"dataset"] == ds))
# 	denom  <- sum(model.data[,"x"]^2)
# 	xt <- model.data[,"x"]/denom
# 	t <- unlist(mclapply(seq_len(R), function(i, model.data, xt, ds.vec) {
#         # model.data[,"x"] <- model.data[,"x"]
#         iny <- sample(model.data[,"y"])
#         # coop::pcor(x,y)
#         y <- standardizeByDatasetInPerm(iny, ds.vec)
#         t0 <- crossprod(xt,y)
# 	}, model.data = model.data, xt = xt, ds.vec = ds.vec, mc.cores=nthread))

#     y <- standardizeByDatasetInPerm(model.data[,"y"], ds.vec)
#     t0 <- crossprod(xt,y)[1]
#     # t0 <- cor(x,y)
#     return(list(t0 = t0, t = t, R = R, pvalue = (sum(abs(t) > abs(t0)) + 1)/(length(t) + 1)))
# }

# current timing, 10 seconds per 1e6, n=320
permFixedEffectC <- function(model.data, R){

	datasets <- unique(model.data[,"dataset"])
	for(ds in datasets){
		myx <- model.data[,"dataset"] == ds
		model.data[myx,"x"] <- scale(model.data[myx,"x"])
	}
	ds.vec <- lapply(datasets, function(ds) return(model.data[,"dataset"] == ds))
	denom  <- sum(model.data[,"x"]^2)
	xt <- model.data[,"x"]/denom
	y <- standardizeByDatasetInPerm(model.data[,"y"], ds.vec)
    t0 <- crossprod(xt,y)[1]

    dsOH <- as.integer(as.numeric(as.factor(model.data[,"dataset"])) - 1)
    DSSize <- as.integer(table(dsOH)[as.character(unique(dsOH))])
    numDS <- as.integer(length(datasets))
    N <- nrow(model.data)
    seed <- as.numeric(runif(2))
    # browser()
    perThread <- ceiling(R/nthread)

    t <- unlist(mclapply(seq_len(nthread), function(i, xt, model.data,dsOH,DSSize,numDS, R, N) {
    	set.seed(seed=NULL)
    	.Call("metaPermC",
                  as.numeric(xt),
                  as.numeric(model.data[,"y"]),
                  as.numeric(t0),
                  dsOH,
                  DSSize,
                  numDS,
                  as.numeric(R),
                  as.numeric(N),
                  as.numeric(runif(2))
                  )}, xt=xt, model.data=model.data,dsOH=dsOH,DSSize=DSSize,numDS=numDS, R=perThread, N=N,mc.cores=nthread))

    return(list("t0" = t0, "t" = t, R = R, pvalue = (sum(abs(t) > abs(t0)) + 1)/(length(t) + 1)))

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

if(!sum(table(model.data[,"dataset"]) > 20) > 2){
    print(paste0("Skipping gene: ", gene))
    
} else {

	corrected_alpha <- corrected_alpha/length(total_gene_list)

	R <- ceiling(1/corrected_alpha)
	## TODO:: Finish below
	# hetTest <- readRDS(file.path(scratch, "perm_meta_hetPermTest", 
	# 	paste0(make.names(drug), "_", make.names(tissue), "_", gene, "_hetPerm_10000_out.rds")))

	# random.effect <- hetTest$pvalue < hetTestCutoff

	# if(random.effect){

	# 	R <- 1e7 

	#     boot.out <- bootMeanEffects(model.data, R)

	#     saveRDS(boot.out, file=file.path(myOutDir, make.names(paste(drug, tissue, gene, "boot", R, "out.RDS", sep="_"))))

	# } else {

		perm.out <- permFixedEffectC(model.data, R)
		perm.out <- c(perm.out, data=model.data)

	    saveRDS(perm.out, file=file.path(myOutDir, make.names(paste("metaPermRes", gene, drug, tissue, "perm.RDS", sep="_"))))

	# }

}




## TODO:
## Here is the approach to take:
## Run bootstrap as follows:
## Sample with replacement studies
## For each study, sample with replacement observations
## For each within study sample, standardize (z-score)
## Fit a random effect model using lme4 (or nlme?) with 0 intercept, random slope per study
## Repeat nboot times
## Use Bootstrap Conf Int, and use the location of effect = 0 to estimate p value
## (This would in effect use a translation null hypothesis - assume the null has similar distribution to observed)
## (Just different mean. )




### forest(meta.res,  xlim=c(-0.2, 0.7), test.overall.fixed=TRUE, test.overall.random=TRUE,comb.fixed=TRUE, comb.random=TRUE,  scientific.pval=FALSE,  layout="meta", calcwidth.hetstat = TRUE, calcwidth.tests=TRUE, rightcols=c("effect", "ci", "w.fixed", "w.random", "pval"), digits.addcols = 6)



