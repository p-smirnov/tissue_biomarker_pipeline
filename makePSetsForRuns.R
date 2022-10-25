library(PharmacoGx)
library(Biobase)
library(data.table)
library(reshape2)
library(CoreGx)
library(SummarizedExperiment)

source('filteringFunctions.R')
## This file does a few things: first of all, it takes unfiltered versoins of psets, and runs
## our qc and uniform concentration range filters on them. 
## It then creates psets for running the pipelines, a set with only rna and a set with only cna
## in them, choosing the prefered profile for each one. these smaller psets are both faster to load
## and allow scripts to just take the only molecular profile in the pset as the one that 
## should be used in a biomarker discovery task. 



## This is set up to run locally, making an assumption that the input psets are in inputDir,
## and the output dir will contain output psets with and rna and cnv folder

inputDir <- "~/Data/unfilteredPSets/"
filteredDir <- "~/Data/TBPInputs/filteredPSets/"

## First, we do filtering



CTRPv2 <- readRDS(file.path(inputDir, "CTRPv2.rds"))

CTRPv2.filtered.sens <- standardizeRawDataConcRange(CTRPv2@sensitivity$info, CTRPv2@sensitivity$raw)

CTRPv2.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=CTRPv2.filtered.sens$sens.raw,
																nthread=20, cap=100, family="normal")
CTRPv2.filtered.profiles <- data.frame("aac_recomputed" = CTRPv2.filtered.profiles.list$AUC, "ic50_recomputed" = CTRPv2.filtered.profiles.list$IC50)
CTRPv2.filtered.profiles.pars <- do.call(rbind,CTRPv2.filtered.profiles.list$pars)
CTRPv2.filtered.profiles.pars <- apply(CTRPv2.filtered.profiles.pars, c(1,2), unlist)

CTRPv2.filtered.profiles <- cbind(CTRPv2.filtered.profiles,CTRPv2.filtered.profiles.pars)

stopifnot(all.equal(rownames(CTRPv2.filtered.profiles), rownames(CTRPv2.filtered.sens$sens.info)))

CTRPv2@sensitivity$info <- CTRPv2.filtered.sens$sens.info
CTRPv2@sensitivity$raw <- CTRPv2.filtered.sens$sens.raw
CTRPv2@sensitivity$profiles <- CTRPv2.filtered.profiles


test <- filterNoisyCurves2(CTRPv2)

CTRPv2@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
CTRPv2@sensitivity$info$noisy.curve <- FALSE
CTRPv2@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE


saveRDS(CTRPv2, file=file.path(filteredDir, "CTRPv2.rds"))




GRAY <- readRDS(file.path(inputDir, "GRAY2017.rds"))

GRAY.filtered.sens <- standardizeRawDataConcRange(GRAY@sensitivity$info, GRAY@sensitivity$raw)

GRAY.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=GRAY.filtered.sens$sens.raw,
																nthread=20, cap=100, family="normal")
GRAY.filtered.profiles <- data.frame("aac_recomputed" = GRAY.filtered.profiles.list$AUC, "ic50_recomputed" = GRAY.filtered.profiles.list$IC50)
GRAY.filtered.profiles.pars <- do.call(rbind,GRAY.filtered.profiles.list$pars)
GRAY.filtered.profiles.pars <- apply(GRAY.filtered.profiles.pars, c(1,2), unlist)

GRAY.filtered.profiles <- cbind(GRAY.filtered.profiles,GRAY.filtered.profiles.pars)

stopifnot(all.equal(rownames(GRAY.filtered.profiles), rownames(GRAY.filtered.sens$sens.info)))

GRAY@sensitivity$info <- GRAY.filtered.sens$sens.info
GRAY@sensitivity$raw <- GRAY.filtered.sens$sens.raw
GRAY@sensitivity$profiles <- GRAY.filtered.profiles



test <- filterNoisyCurves2(GRAY)

GRAY@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
GRAY@sensitivity$info$noisy.curve <- FALSE
GRAY@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(GRAY, file=file.path(filteredDir, "GRAY2017.rds"))





GDSC1 <- readRDS(file.path(inputDir, "GDSC1.rds"))

GDSC1.filtered.sens <- standardizeRawDataConcRange(GDSC1@sensitivity$info, GDSC1@sensitivity$raw)

GDSC1.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=GDSC1.filtered.sens$sens.raw,
																nthread=12, cap=100, family="normal")

GDSC1.filtered.profiles <- data.frame("aac_recomputed" = GDSC1.filtered.profiles.list$AUC, "ic50_recomputed" = GDSC1.filtered.profiles.list$IC50)
GDSC1.filtered.profiles.pars <- do.call(rbind,GDSC1.filtered.profiles.list$pars)
GDSC1.filtered.profiles.pars <- apply(GDSC1.filtered.profiles.pars, c(1,2), unlist)

GDSC1.filtered.profiles <- cbind(GDSC1.filtered.profiles,GDSC1.filtered.profiles.pars)

stopifnot(all.equal(rownames(GDSC1.filtered.profiles), rownames(GDSC1.filtered.sens$sens.info)))

GDSC1@sensitivity$info <- GDSC1.filtered.sens$sens.info
GDSC1@sensitivity$raw <- GDSC1.filtered.sens$sens.raw
GDSC1@sensitivity$profiles <- GDSC1.filtered.profiles



test <- filterNoisyCurves2(GDSC1)

GDSC1@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
GDSC1@sensitivity$info$noisy.curve <- FALSE
GDSC1@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(GDSC1, file=file.path(filteredDir, "GDSC1.rds"))




## 



GDSC2 <- readRDS(file.path(inputDir, "GDSC2.rds"))

GDSC2.filtered.sens <- standardizeRawDataConcRange(GDSC2@sensitivity$info, GDSC2@sensitivity$raw)

GDSC2.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=GDSC2.filtered.sens$sens.raw,
																nthread=20, cap=100, family="normal")

GDSC2.filtered.profiles <- data.frame("aac_recomputed" = GDSC2.filtered.profiles.list$AUC, "ic50_recomputed" = GDSC2.filtered.profiles.list$IC50)
GDSC2.filtered.profiles.pars <- do.call(rbind,GDSC2.filtered.profiles.list$pars)
GDSC2.filtered.profiles.pars <- apply(GDSC2.filtered.profiles.pars, c(1,2), unlist)

GDSC2.filtered.profiles <- cbind(GDSC2.filtered.profiles,GDSC2.filtered.profiles.pars)

stopifnot(all.equal(rownames(GDSC2.filtered.profiles), rownames(GDSC2.filtered.sens$sens.info)))

GDSC2@sensitivity$info <- GDSC2.filtered.sens$sens.info
GDSC2@sensitivity$raw <- GDSC2.filtered.sens$sens.raw
GDSC2@sensitivity$profiles <- GDSC2.filtered.profiles



test <- filterNoisyCurves2(GDSC2)

GDSC2@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
GDSC2@sensitivity$info$noisy.curve <- FALSE
GDSC2@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(GDSC2, file=file.path(filteredDir, "GDSC2.rds"))





## 



gCSI <- readRDS(file.path(inputDir, "gCSI2.rds"))

gCSI.filtered.sens <- standardizeRawDataConcRange(gCSI@sensitivity$info, gCSI@sensitivity$raw)

gCSI.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=gCSI.filtered.sens$sens.raw,
																nthread=12, cap=100, family="normal")

gCSI.filtered.profiles <- data.frame("aac_recomputed" = gCSI.filtered.profiles.list$AUC, "ic50_recomputed" = gCSI.filtered.profiles.list$IC50)
gCSI.filtered.profiles.pars <- do.call(rbind,gCSI.filtered.profiles.list$pars)
gCSI.filtered.profiles.pars <- apply(gCSI.filtered.profiles.pars, c(1,2), unlist)

gCSI.filtered.profiles <- cbind(gCSI.filtered.profiles,gCSI.filtered.profiles.pars)

stopifnot(all.equal(rownames(gCSI.filtered.profiles), rownames(gCSI.filtered.sens$sens.info)))

gCSI@sensitivity$info <- gCSI.filtered.sens$sens.info
gCSI@sensitivity$raw <- gCSI.filtered.sens$sens.raw
gCSI@sensitivity$profiles <- gCSI.filtered.profiles



test <- filterNoisyCurves2(gCSI)

gCSI@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
gCSI@sensitivity$info$noisy.curve <- FALSE
gCSI@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(gCSI, file=file.path(filteredDir, "gCSI2018.rds"))





CCLE <- readRDS(file.path(inputDir, "CCLE.rds"))

CCLE.filtered.sens <- standardizeRawDataConcRange(CCLE@sensitivity$info, CCLE@sensitivity$raw)

CCLE.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=CCLE.filtered.sens$sens.raw,
																nthread=12, cap=100, family="normal")

CCLE.filtered.profiles <- data.frame("aac_recomputed" = CCLE.filtered.profiles.list$AUC, "ic50_recomputed" = CCLE.filtered.profiles.list$IC50)
CCLE.filtered.profiles.pars <- do.call(rbind,CCLE.filtered.profiles.list$pars)
CCLE.filtered.profiles.pars <- apply(CCLE.filtered.profiles.pars, c(1,2), unlist)

CCLE.filtered.profiles <- cbind(CCLE.filtered.profiles,CCLE.filtered.profiles.pars)

stopifnot(all.equal(rownames(CCLE.filtered.profiles), rownames(CCLE.filtered.sens$sens.info)))

CCLE@sensitivity$info <- CCLE.filtered.sens$sens.info
CCLE@sensitivity$raw <- CCLE.filtered.sens$sens.raw
CCLE@sensitivity$profiles <- CCLE.filtered.profiles



test <- filterNoisyCurves2(CCLE)

CCLE@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
CCLE@sensitivity$info$noisy.curve <- FALSE
CCLE@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(CCLE, file=file.path(filteredDir, "CCLE.rds"))





UHNBreast <- readRDS(file.path(inputDir, "UHNBreast.rds"))

UHNBreast.filtered.sens <- standardizeRawDataConcRange(UHNBreast@sensitivity$info, UHNBreast@sensitivity$raw)

UHNBreast.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=UHNBreast.filtered.sens$sens.raw,
																nthread=12, cap=100, family="normal")

UHNBreast.filtered.profiles <- data.frame("aac_recomputed" = UHNBreast.filtered.profiles.list$AUC, "ic50_recomputed" = UHNBreast.filtered.profiles.list$IC50)
UHNBreast.filtered.profiles.pars <- do.call(rbind,UHNBreast.filtered.profiles.list$pars)
UHNBreast.filtered.profiles.pars <- apply(UHNBreast.filtered.profiles.pars, c(1,2), unlist)

UHNBreast.filtered.profiles <- cbind(UHNBreast.filtered.profiles,UHNBreast.filtered.profiles.pars)

stopifnot(all.equal(rownames(UHNBreast.filtered.profiles), rownames(UHNBreast.filtered.sens$sens.info)))

UHNBreast@sensitivity$info <- UHNBreast.filtered.sens$sens.info
UHNBreast@sensitivity$raw <- UHNBreast.filtered.sens$sens.raw
UHNBreast@sensitivity$profiles <- UHNBreast.filtered.profiles



test <- filterNoisyCurves2(UHNBreast)

UHNBreast@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
UHNBreast@sensitivity$info$noisy.curve <- FALSE
UHNBreast@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(UHNBreast, file=file.path(filteredDir, "UHNBreast.rds"))





PRISM <- readRDS(file.path(inputDir, "PRISM.rds"))

PRISM.filtered.sens <- standardizeRawDataConcRange(PRISM@sensitivity$info, PRISM@sensitivity$raw)

PRISM.filtered.profiles.list <- PharmacoGx:::.calculateFromRaw(raw.sensitivity=PRISM.filtered.sens$sens.raw,
																nthread=12, cap=100, family="normal")

PRISM.filtered.profiles <- data.frame("aac_recomputed" = PRISM.filtered.profiles.list$AUC, "ic50_recomputed" = PRISM.filtered.profiles.list$IC50)
PRISM.filtered.profiles.pars <- do.call(rbind,PRISM.filtered.profiles.list$pars)
PRISM.filtered.profiles.pars <- apply(PRISM.filtered.profiles.pars, c(1,2), unlist)

PRISM.filtered.profiles <- cbind(PRISM.filtered.profiles,PRISM.filtered.profiles.pars)

stopifnot(all.equal(rownames(PRISM.filtered.profiles), rownames(PRISM.filtered.sens$sens.info)))

PRISM@sensitivity$info <- PRISM.filtered.sens$sens.info
PRISM@sensitivity$raw <- PRISM.filtered.sens$sens.raw
PRISM@sensitivity$profiles <- PRISM.filtered.profiles



test <- filterNoisyCurves2(PRISM)

PRISM@sensitivity$profiles[test$noisy,c("aac_recomputed","ic50_recomputed")] <- NA_real_
PRISM@sensitivity$info$noisy.curve <- FALSE
PRISM@sensitivity$info[test$noisy,"noisy.curve"] <- TRUE



saveRDS(PRISM, file=file.path(filteredDir, "PRISM.rds"))

