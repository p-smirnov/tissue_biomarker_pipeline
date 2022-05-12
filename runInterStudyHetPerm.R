library(PharmacoGx)
library(meta)
library(SummarizedExperiment)
library(coop)
library(foreach)
# library(doParallel)
# library(doRNG)
library(iterators)
library(lme4)
library(data.table)

library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
RhpcBLASctl::blas_set_num_threads(1)
#registerDoParallel(40)

args <- commandArgs(trailingOnly = TRUE)

method <- "perm"


drug <- args[1]

# tissue <- "Breast"
tissue <- args[2]

gene <- args[3]

badchars <- "[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[ ]|[(]|[)]"

make.names.2 <- function(x) return(gsub(pat=badchars, rep=".", x))



R  <- as.numeric(args[4]) ## TODO: need method to pick this


home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")
dataDir <- Sys.getenv("DATA")

myDataDir <- dataDir
mySigDir <- file.path(project, paste0("pearson_",method,"_res"))
myOutDir <- file.path(project, paste0("hetTest"))
runlistDir <- file.path(project, paste0("runlist_files"))

print(drug)
print(tissue)
print(gene)

containername <- Sys.getenv("containername", unset=NA_character_)
snakemake <- as.numeric(Sys.getenv("SNAKEMAKE", unset=0))


if(!is.na(containername)){
	myOutDir <- file.path(containername, myOutDir)
	mySigDir <- file.path(containername, mySigDir)
	myDataDir <- file.path(containername, myDataDir)
	project <- file.path(containername, project)
	runlistDir <- file.path(containername, runlistDir)


}


loadPSet <- function(psetName, tissue){


	switch(psetName, 
		   CCLE = {
		   		pset <- readRDS(file.path(myDataDir,"CCLE.rds"))
		   	}, CCLE.CTRPv2 = {
		   		pset <- readRDS(file.path(myDataDir,"CCLE.CTRPv2.rds"))
		   	}, CCLE.PRISM = {
	   			pset <- readRDS(file.path(myDataDir,"CCLE.PRISM.rds"))
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



toRunExtended <- data.frame(fread(file.path(runlistDir,"toRunMetaByGene.txt"), header=FALSE))

# need to do this "trick" because names are made path safe, and arguments are derived from paths for snakemake's sake 

if(snakemake){
  drug <- unique(toRunExtended[,3])[make.names.2(unique(toRunExtended[,3])) == drug]
  tissue <- unique(toRunExtended[,2])[make.names.2(unique(toRunExtended[,2])) == tissue]
  gene <- unique(toRunExtended[,1])[make.names.2(unique(toRunExtended[,1])) == gene]
}



toRunThis <- toRunExtended[toRunExtended[,3] == drug & toRunExtended[,2] == tissue & toRunExtended[,1] == gene, ]
# drug <- unique(toRunThis[,3])
# tissue <- unique(toRunThis[,2])
# gene <- unique(toRunThis[,1])
pSets <- toRunThis[,4]



if(!file.exists(myOutDir)) dir.create(myOutDir)


pset.list <- lapply(pSets, loadPSet, tissue=tissue)
names(pset.list) <- pSets

runPermutations <- function(beta_obs, errors, model.data, R){

	perm_res <- numeric(R)

	perm_res <- foreach(i =icount(R), .combine=c, .inorder=FALSE) %do% {

		errors_p <- sample(errors)

		model.data$y_p <- beta_obs*model.data[,"x"] + errors_p

		lmer_p <- suppressMessages(lmer(y_p~(x + 0| dataset) + x + 0, model.data))

		sum(unlist(ranef(lmer_p))^2)

	}
	return(perm_res)
}


psetMolProf <- lapply(pset.list, function(pset.pars){
	
	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	mData <- pset.pars$mData

	mol.prof <- assay(summarizeMolecularProfiles(pset, mData, cell.lines = chosen.cells))
	return(mol.prof)
})

names(psetMolProf) <- pSets


y.list <- lapply(pset.list, function(pset.pars){

	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	
	drug.res <- summarizeSensitivityProfiles(pset, "aac_recomputed", cell.lines = chosen.cells)[drug,]


	return(drug.res)


})


x.list <- lapply(pset.list, function(pset.pars){

	pset <- pset.pars$pset
	chosen.cells <- pset.pars$chosen.cells
	mData <- pset.pars$mData

	mol.prof <- psetMolProf[[name(pset)]]

	if(grepl(pat="ENSG", x=gene)){
		myx <- grep(paste0("^", gene, "(\\.[0-9]+)?$"), rownames(mol.prof))
	} else {
		myx <- grep(paste0("^", gene, "$"), rownames(mol.prof))
	}

	if(!length(myx)) return(rep(NA_real_, times=length(pset.pars$chosen.cells)))
	if(length(myx) > 1) stop(paste0("multiple genes matched:", gene, " id in pset ", name(pset), ", please investigate"))

	mol.prof <- mol.prof[myx,]

	return(mol.prof)


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


model.data <- regressOutTissueAndScale(model.data)

if(sum(table(model.data[,"dataset"]) >= 20) > 2){
	lmer1 <- lmer(y~(x + 0| dataset) + x + 0, model.data)
	tobs <- sum(unlist(ranef(lmer1))^2)

	beta_obs <- fixef(lmer1)[1]

	errors <- model.data[,"y"] - fixef(lmer1)[1]*model.data[,"x"]
	permRes <- runPermutations(beta_obs, errors, model.data, R)
	heterogeneityRes <- list(gene = gene, pvalue = mean(permRes > tobs))




	saveRDS(heterogeneityRes, file=file.path(myOutDir, make.names.2(paste(drug, tissue, gene, "hetPerm", R, "out.rds", sep="_"))))

} else {
	message(paste0("Skipping gene: ", gene))
}



# Sirolimus_Myeloid_ENSG00000164850
# 5.Fluorouracil_Myeloid_ENSG00000112697
# Navitoclax_Lung_ENSG00000002586

# Lapatinib_Breast_ENSG00000197520
# ENSG00000100300 Breast Alpelisib CCLE.CTRPv2
