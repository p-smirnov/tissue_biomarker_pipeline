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
library(boot)



library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
RhpcBLASctl::blas_set_num_threads(1)
# registerDoParallel(40)



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

myDataDir <- Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myOutDir <- file.path(project, paste0(method, "_meta_out"))
myPvalDir <- file.path(project, paste0(method, "_meta_perm_sig"))


runDir <- args[6]

containername <- Sys.getenv("containername", unset=NA_character_)

if(!is.na(containername)){
    myDataDir <- file.path(containername, myDataDir)
    mySigDir <- file.path(containername, mySigDir)
    myOutDir <- file.path(containername, myOutDir)
    myPvalDir <- file.path(containername, myPvalDir)
    runDir <- file.path(containername, runDir)


}






codeDir <- args[7]
dyn.load(file.path(codeDir,"metaPermC.so"))
dyn.load(file.path(codeDir,"metaPermCTissue.so"))


badchars <- "[ ]|[/]|[:]|[-]"

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
		}, Tavor = {
			pset <- readRDS(file.path(myDataDir, "Tavor.rds"))
		}, BeatAML = {
			pset <- readRDS(file.path(myDataDir, "BeatAML.rds"))
		}, "FIMM-AML-MCM" = {
			pset <- readRDS(file.path(myDataDir, "FIMM_MCM.rds"))
		}, {stop("Please Provide a valid pset")})

	mData <- mDataNames(pset)


	gene_type_col <- ifelse("GeneBioType" %in% colnames(featureInfo(pset, mData)), "GeneBioType", "gene_type")

	ft <- rownames(featureInfo(pset, mData))[which(featureInfo(pset, mData)[[gene_type_col]] == "protein_coding")]

	if(is.na(tissue)||tissue=="all"){
		chosen.cells <- cellNames(pset)
	} else {
		chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
	}

	return(list(pset = pset, mData = mData, chosen.cells = chosen.cells))
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



if(!file.exists(myOutDir)) dir.create(myOutDir)
if(!file.exists(myPvalDir)) dir.create(myPvalDir)


pset.list <- lapply(pSets, loadPSet, tissue=tissue)
names(pset.list) <- pSets


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
# 			model.data[myx,"x"] <- scale(model.data[myx,"x"],scale=T)
# 			model.data[myx,"y"] <- scale(model.data[myx,"y"],scale=T)	

# 		}
# 		return(cor(model.data[,"x"], model.data[,"y"]))
# 	}, model.data = model.data, datasets = datasets, mc.cores=nthread))
# 		for(ds in datasets){
# 			myx <- model.data[,"dataset"] == ds
# 			model.data[myx,"x"] <- scale(model.data[myx,"x"],scale=T)
# 			model.data[myx,"y"] <- scale(model.data[myx,"y"],scale=T)	

# 		}
#      t0 <- NA_real_
#     t0 <- cor(model.data[,"x"],model.data[,"y"])
#     return(list(t0 = t0, t = t, R = R, pvalue = (sum(abs(t) > abs(t0)) + 1)/(length(t) + 1)))
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
# 	denom  <- sum(model.data[,"x"]^2) ### Note that this just ends up calculating N-M, where N is num samples, M is num groups
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
##FIXME:: not properly accounting for equality up to numerical precision here!
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
    # DSSize <- as.integer(table(dsOH)[as.character(unique(dsOH))])
    DSSize <- as.integer(table(dsOH))

    numDS <- as.integer(length(datasets))
    N <- nrow(model.data)
    seed <- as.numeric(runif(2))
    # browser()
    perThread <- ceiling(R/nthread)

    t <- unlist(mclapply(seq_len(nthread), function(i, xt, t0, model.data,dsOH,DSSize,numDS, R, N) {
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
                  )}, xt=xt, t0=t0, model.data=model.data,dsOH=dsOH,DSSize=DSSize,numDS=numDS, R=perThread, N=N,mc.cores=nthread))

    return(list("t0" = t0, "t" = t, R = R, pvalue = (sum(abs(t) >= abs(t0)) + 1)/(length(t) + 1)))

}

# dyn.load("metaPermCTissue.so")
# dyn.load("metaPermC.so")
# library(parallel)


# nthread=1
# model.data <- data.frame(x=runif(1000), y=runif(1000), dataset=as.character(sample(c("A","B","C","D","E"), 1000, replace=TRUE)),
#  tissueid=as.character(sample(c("a", "b", "c", "d", "e", "f","g","h","i"), 1000, replace=TRUE)))

# permFixedEffectC(model.data, 100)

standardizeByDatasetInPerm <- function(iny,ds.vec){
	# y <- numeric(length(iny))
	for(myx in ds.vec){
		# model.data[myx,"x"] <- scale(model.data[myx,"x"])
		iny[myx] <- scale(iny[myx])
	}
	return(iny)
}


regressOutTissueAndScale <- function(model.data){

	for(ds in model.data[,"dataset"]){
		myx <- model.data[,"dataset"]==ds
		if(length(unique(model.data[myx,"tissueid"]))>=2){

			model.data[myx,"x"] <- residuals(lm(x~tissueid, data=model.data[myx,]))
			model.data[myx,"y"] <- residuals(lm(y~tissueid, data=model.data[myx,]))


		}
		model.data[myx,"x"] <- scale(model.data[myx, "x"])
		model.data[myx,"y"] <- scale(model.data[myx, "y"])

	}
	return(model.data)
}


permFixedEffectCTissue <- function(model.data, R){

	
		datasets <- unique(model.data[,"dataset"])
		tissues <- unique(model.data[,"tissueid"])

    dsOH <- as.integer(as.numeric(as.factor(model.data[,"dataset"])) - 1)
    DSSize <- as.integer(table(dsOH))
    # browser()
    ## TODO: make sure this TissueSize works properly
    ## This does not work properly!!! Need to fix!!!
    tissueOH <- as.integer(as.numeric(as.factor(model.data[,"tissueid"]))-1)
    # TissueSize <- as.integer(table(tissueOH))
    TissueSize <- as.vector(table(tissueOH, dsOH))
    # browser()
    ## TODO: FIXME: Need to make sure this is done per dataset IN DATASET ORDER!
    ## should be size numDS*numTissue 

    numDS <- as.integer(length(datasets))
		numTissue <- as.integer(length(tissues))



		# ## Calculating t0
		# res.x <- residuals(lm(x~tissueid, data=model.data))
		# res.y <- residuals(lm(y~tissueid, data=model.data))
		# ds.vec <- lapply(datasets, function(ds) return(model.data[,"dataset"] == ds))
		# # print(res.y)

		# # browser()

		# yt <- standardizeByDatasetInPerm(res.y, ds.vec)
		# xt <- standardizeByDatasetInPerm(res.x, ds.vec)
		# browser()
		model.data.mod <- regressOutTissueAndScale(model.data)
		xt <- model.data.mod[,"x"]
		yt <- model.data.mod[,"y"]
		t0 <- cor(xt, yt)

    N <- nrow(model.data)
    seed <- as.numeric(runif(2))
    # browser()
    perThread <- ceiling(R/nthread)

    t <- unlist(mclapply(seq_len(nthread), function(i, 
    																								model.data, 
    																								dsOH,
    																								DSSize,
    																								numDS,
    																								tissueOH,
    																								TissueSize,
    																								numTissue,
    																								t0,
    																								R, 
    																								N) {
    	set.seed(seed=NULL)
    	.Call("metaPermCTissue",
                  as.numeric(model.data[,"x"]),
                  as.numeric(model.data[,"y"]),
                  as.numeric(t0),
                  dsOH,
                  DSSize,
                  tissueOH,
                  TissueSize,
                  numDS,
                  numTissue,
                  as.numeric(R),
                  as.numeric(N),
                  as.numeric(runif(2))
                  )
    },  model.data=model.data, dsOH=dsOH, DSSize=DSSize, numDS=numDS, tissueOH=tissueOH, TissueSize=TissueSize, 
    		numTissue=numTissue, t0=t0, R=perThread, N=N, mc.cores=nthread))
    return(list("t0" = t0, "t" = t, R = R, pvalue = (sum(abs(t) >= abs(t0)) + 1)/(length(t) + 1)))

}

# system.time(permFixedEffectCTissue(model.data, 1000000))

x.list <- lapply(pset.list, function(pset.pars){

	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	mData <- pset.pars$mData

	mol.prof <- assay(summarizeMolecularProfiles(pset, mData, cell.lines = chosen.cells))

	if(grepl(pat="ENSG", x=gene)){
			myx <- grep(paste0("^", gene, "(\\.[0-9]+)?$"), rownames(mol.prof))
	} else {
			myx <- grep(paste0("^", gene, "$"), rownames(mol.prof))
	}	
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

tissue.list <- lapply(pset.list, function(pset.pars){

	return(cellInfo(pset.pars$pset)[pset.pars$chosen.cells,"tissueid"])

})

model.data <- data.frame(x = unlist(x.list),
						 y = unlist(y.list),
						 dataset = unlist(lapply(pset.list, function(pset.pars){
						 	return(rep(name(pset.pars$pset), times=length(pset.pars$chosen.cells)))
						 	})),
						 tissueid = unlist(tissue.list)
						)


model.data <- model.data[complete.cases(model.data),]


if(!sum(table(model.data[,"dataset"]) > 20) > 2){
    print(paste0("Skipping gene: ", gene))
    
} else {

	corrected_alpha <- corrected_alpha/length(total_gene_list)

	R <- ceiling(1/corrected_alpha)


	if(length(unique(model.data$tissueid))==1){
		perm.out <- permFixedEffectC(model.data, R)
		perm.out <- c(perm.out, data=model.data)		
	} else {
		perm.out <- permFixedEffectCTissue(model.data, R)
		perm.out <- c(perm.out, data=model.data)
	}

	saveRDS(perm.out, file=file.path(myOutDir, make.names(paste("metaPermRes", gene, drug, tissue, "perm.RDS", sep="_"))))

}


## Doing bootstrapping for Conf Int estimation


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

## techinically above is a subset, but branching is efficient
getBootSampleTissue <- function(model.data, scale=TRUE){

  sampled.datasets <- sample(unique(model.data$dataset), replace=TRUE)
  
  names(sampled.datasets) <- make.unique(sampled.datasets)

  sampled.data <- lapply(names(sampled.datasets), function(dt){

    sm.dt <- model.data[sample(which(model.data$dataset == sampled.datasets[dt]), replace=TRUE), ]

    if(length(unique(sm.dt$tissueid))>1){
    	x <- as.vector(residuals(lm(x~tissueid, data=sm.dt))) ## Otherwise they get dimnames???
    	y <- as.vector(residuals(lm(y~tissueid, data=sm.dt)))
    } else {
    	x <- sm.dt[,"x"]
    	y <- sm.dt[,"y"]
    }
    sm.dt$dataset <- dt
        sm.dt[,1] <- scale(x)
        sm.dt[,2] <- scale(y)
    return(sm.dt)
    })

  sampled.data <- do.call(rbind, sampled.data)

  return(sampled.data)
}

bootCor <- function(model.data,R, t0=perm.out$t0){

    N <- nrow(model.data)

    t <- unlist(mclapply(seq_len(R), function(i, model.data) {

    	if(length(unique(model.data$tissueid))==1){
    		model.data.sampled <- getBootSample(model.data)
    	} else {
    		model.data.sampled <- getBootSampleTissue(model.data)
    	}
      denom  <- sum(model.data.sampled[,"x"]^2)
      xt <- model.data.sampled[,"x"]/denom


      crossprod(xt,model.data.sampled[,"y"])[1]
    }, model.data = model.data, mc.cores=nthread))

    #  model.data.scaled <- do.call(rbind, lapply(model.data$dataset, function(dt){
    #     sm.dt <- model.data[model.data$dataset == dt,]
    #     sm.dt[,1] <- scale(sm.dt[,1])
    #     sm.dt[,2] <- scale(sm.dt[,2])
    #     return(sm.dt)
    #     }))

    # denom  <- sum(model.data.scaled[,"x"]^2)
    # xt <- model.data.scaled[,"x"]/denom

    # t0 <- crossprod(xt, model.data.scaled[,"y"])[1]


  dim(t) <- c(R,1)
  res <- list(
      t0 = t0,
      t = t,
      R = R,
      data = model.data,
      seed = NA_real_,
      sim = "ordinary",
      stype = "i",
      call = match.call(),
      strata = rep(1, N),
      weights = rep(1/N, N)
  )
    attr(res, "class") <- "boot"
    attr(res, "boot_type") <- "boot"
    return(res)
}



# perm.out <- readRDS(myInFl)

# p.value <- mean(abs(perm.out$t0) <= abs(perm.out$t))
# if(p.value == 0){
#   p.value <- 1/(length(perm.out$t)+1)
# }

p.value <- perm.out$pvalue

# model.data <- data.frame(x=perm.out$data.x, y=perm.out$data.y,dataset=perm.out$data.dataset)

boot.res <- bootCor(model.data, 10000) 

ci.95.boot <- boot.ci(boot.res, conf = 0.95, type="perc")

perm.sig.res <- list("p.value" = p.value, ci = ci.95.boot, LR = R, R = length(perm.out$t))

saveRDS(perm.sig.res, file=file.path(myPvalDir, make.names(paste("permSig", gene, drug, tissue, "out.rds", sep="_"))))


