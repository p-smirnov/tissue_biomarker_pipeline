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



args <- commandArgs(trailingOnly = TRUE)
# myToRunFileName <- args[[1]]

method <- "perm"
hetTestCutoff <- 0.1
## ACSS2 is interesting



# psetName <- args[1]

# drug <- "Lapatinib"

## Asking for 100 times more permutations than the cutoff for alpha 
corrected_alpha <- 0.05/100

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")


dataDir <- Sys.getenv("DATA")

myDataDir <- dataDir


myOutDir <- args[[2]]
inDir <- args[[1]]




containername <- Sys.getenv("containername", unset=NA_character_)

if(!is.na(containername)){
    myOutDir <- file.path(containername, myOutDir)
    inDir <- file.path(containername, inDir)
    myDataDir <- file.path(containername, myDataDir)


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

    # if(is.na(tissue)){
    #   chosen.cells <- cellNames(pset)
    # } else {
    #   chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
    # }

    return(list(pset = pset, mData = mData))
}


toRunMetaRes <- fread(file=file.path(inDir, "metaHetTestRes.txt"))

toRunByGene <- toRunMetaRes[hetTestRes == TRUE]


myToRunFileName <- file.path(inDir, "toRunMetaByGene.txt")

toRunExtended <- fread(myToRunFileName, header=FALSE)

colnames(toRunExtended) <- c("Gene", "Tissue", "Drug", "PSet", "Significant")

toRunMetaRes <- merge(toRunByGene, toRunExtended, on=c("Gene", "Tissue", "Drug"))

total_gene_list <- toRunExtended[,unique(Gene)]


if(!file.exists(myOutDir)) dir.create(myOutDir, recursive = TRUE)

pSets <- unique(toRunExtended[,unique(PSet)])

pset.list <- lapply(pSets, loadPSet)
names(pset.list) <- pSets

# toRunByGene <- read.csv("~/toRunMetaByGeneBoot.txt", header=FALSE)

mol.list <- lapply(pset.list, function(pset.pars){

    pset <- pset.pars$pset
    mData <- pset.pars$mData

    mol.prof <- assay(summarizeMolecularProfiles(pset, mData))

    return(mol.prof)

})

names(mol.list) <- names(pset.list)
sens.list <- lapply(pset.list, function(pset.pars){

    pset <- pset.pars$pset
        
    drug.res <- summarizeSensitivityProfiles(pset, "aac_recomputed")

    return(drug.res)
})

names(sens.list) <- names(pset.list)



# gene.fls <- list.files(path="~/featureLists/", full.names=TRUE)
# gene.fls.short <- list.files(path="~/featureLists/")

# gene.fls.list <- lapply(gene.fls, readLines)
# names(gene.fls.list) <- gene.fls.short
# #1:nrow(toRunByGene)

## TODO:: check that I don't need to filter to only valid tissues here!
## TODO:: test me
for(i in seq_len(nrow(toRunByGene))){

    drug <- toRunByGene[i,Drug]
    tissue <- toRunByGene[i,Tissue]
    gene <- toRunByGene[i,Gene]



    pSetsToRun <- toRunMetaRes[Drug == drug & Tissue == tissue & Gene == gene, PSet]
    pset.listLocal <- pset.list[pSetsToRun]


    x.list <- lapply(pset.listLocal, function(pset.pars){

        pset <- pset.pars$pset
        if(is.na(tissue)||tissue=="all"){
            chosen.cells <- cellNames(pset)
        } else{
            chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
        }
        mData <- pset.pars$mData

        mol.prof <- mol.list[[name(pset)]]
        mol.prof <- mol.prof[,chosen.cells]

        if(grepl(pat="ENSG", x=gene)){
            myx <- grep(paste0("^", gene, "(\\.[0-9]+)?$"), rownames(mol.prof))
        } else {
            myx <- grep(paste0("^", gene, "$"), rownames(mol.prof))
        }  

        if(!length(myx)) return(rep(NA_real_, times=length(chosen.cells)))
        if(length(myx) > 1) stop(paste0("multiple genes matched:", gene, " id. please investigate"))

        mol.prof <- mol.prof[myx,]

        return(mol.prof)


        })

    y.list <- lapply(pset.listLocal, function(pset.pars){

        pset <- pset.pars$pset
        if(is.na(tissue)||tissue=="all"){
            chosen.cells <- cellNames(pset)
        } else{
            chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
        }        
        drug.res <- sens.list[[name(pset)]][drug,chosen.cells]


        return(drug.res)


        })

    tissue.list <- lapply(pset.listLocal, function(pset.pars){

        pset <- pset.pars$pset
        if(is.na(tissue)||tissue=="all"){
            tissue <- cellInfo(pset)$tissueid
        } else{
            tissue <- cellInfo(pset)[which(cellInfo(pset)$tissueid == tissue),"tissueid"]
        }        

        return(tissue)



        })



    dataset.list <- (lapply(pset.listLocal, function(pset.pars){
        pset <- pset.pars$pset
        if(is.na(tissue)||tissue=="all"){
            chosen.cells <- cellNames(pset)
        } else{
            chosen.cells <- cellNames(pset)[which(cellInfo(pset)$tissueid == tissue)]
        }         
        return(rep(name(pset.pars$pset), times=length(chosen.cells)))
    
    }))

    model.data <- data.frame(x = unlist(x.list),
                             y = unlist(y.list),
                             dataset = unlist(dataset.list),
                             tissueid = unlist(tissue.list))


    model.data <- model.data[complete.cases(model.data),]


    # gene.fls.this <- paste0(make.names(tissue), ".", make.names(pSetsToRun), ".csv")
    # total_gene_list <- .unionList(gene.fls.list[gene.fls.this])

    corrected_alpha.this <- corrected_alpha/length(total_gene_list)

    R <- ceiling(1/corrected_alpha.this)
    model.data$R <- R


    write.csv(model.data, file=file.path(myOutDir,  make.names(paste0("modelData_", gene, "_", drug, "_", tissue, ".csv"))))
    if(!i%%100) print(i)
}






