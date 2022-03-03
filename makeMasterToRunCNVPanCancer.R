## This file operates on assumption that if the drug-cnv gene combo existed in the per-tissue 
## analysis, it should also be investigated for the pan-cancer analysis


library(PharmacoGx)
library(data.table)

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(home, "Data", "TBPInputs", "cnv")

all.dt.pertissue <- fread(file.path(scratch, "cnvMasterToRunList.txt"), header=F)

all.dt.pertissue[,V2:="all"]

all.dt <- unique(all.dt.pertissue)


write.table(all.dt, file.path(scratch, "cnvPancancerMasterToRunList.txt"), quote=FALSE, row.names=FALSE, sep=",", col.names=FALSE)



