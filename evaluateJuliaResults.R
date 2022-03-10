library(PharmacoGx)
library(meta)
library(SummarizedExperiment)
library(coop)
# library(foreach)
# library(doParallel)
# library(doRNG)
library(iterators)
library(lme4)
library(boot)
library(data.table)


library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
# RhpcBLASctl::blas_set_num_threads(1)
# registerDoParallel(40)

bootR <- "10000000"

args <- commandArgs(trailingOnly = TRUE)

method <- "perm"
hetTestCutoff <- 0.1
## ACSS2 is interesting



# psetName <- args[1]

# drug <- "Lapatinib"

drug <- args[1]

# tissue <- "Breast"
tissue <- args[2]

gene <- args[3]


## 206 drugs/tissues being tested + 2 orders of magnitude for accurate computations
R <- as.numeric(args[4]) ## In this case, its the number of samples to use for L estimation

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")

# myBootDir <- file.path(scratch, paste0(method, "_meta_bootstrap"))
myPvalDir <- file.path(project, paste0(method, "_meta_boot_sig"))

myInFl <- args[5]
myDataFl <- args[6]

if(!file.exists(myPvalDir)) dir.create(myPvalDir)

myBootDir <- file.path(scratch, "juliaBoot")


containername <- Sys.getenv("containername", unset=NA_character_)

if(!is.na(containername)){
    myBootDir <- file.path(containername, myBootDir)
    myPvalDir <- file.path(containername, myPvalDir)
    # myDataFl <- file.path(containername, myDataFl)
    # myInFl <- file.path(containername, myInFl)


}



model.data <- fread(myDataFl)

model.data[,R:=NULL]
model.data[,V1:=NULL]

### First we handle the CI and p-value computation for the bootstrapped meta-estimate

parseBootResults <- function(unparsed){

  R <- as.numeric(unparsed[[match("R:", unparsed) + 1]])
  Tstart <- match("t:", unparsed) + 1
  Tend <- Tstart + R - 1
  t0 <- as.numeric(unparsed[[2]]) ## NB:: these values are incorrect!!!! Code has been fixed, but rerunning is costly!
  t <- as.numeric(unlist(unparsed[Tstart:Tend]))
  N <- as.numeric(unparsed[[4]])
  dim(t) <- c(R,1)
  res <- list(
      t0 = t0,
      t = t,
      R = R,
      data = as.data.frame(model.data),
      seed = NA_real_,
      sim = "ordinary",
      stype = "i",
      call = match.call(),
      strata = NA,
      weights = rep(1/N, N)
  )
  attr(res, "class") <- "boot"
  attr(res, "boot_type") <- "boot"
  return(res)

}


# randomEffectBootAndSampleArray <- function(model.data, R){

# 	SampleArrayList <- list()

# 	for(dataset in unique(model.data$dataset)){

# 		SampleArrayList[[dataset]] <- matrix(0,nrow=R, ncol=sum(model.data$dataset==dataset))
# 	}

# 	t <- numeric(R)

# 	for(i in seq_len(R)){

# 		sampled.datasets <- sample(unique(model.data$dataset), replace=TRUE)
	
# 		names(sampled.datasets) <- make.unique(sampled.datasets)

# 		sampled.data <- list()
# 		for(dt in names(sampled.datasets)) {

# 			myDtX <- model.data$dataset == sampled.datasets[dt]
# 			model.data.dt <- model.data[myDtX, ]
# 			myx <- sample.int(nrow(model.data.dt), replace=TRUE)
# 			SampleArrayList[[sampled.datasets[dt]]][i,] <- tabulate(myx, ncol(SampleArrayList[[sampled.datasets[dt]]])) 

# 			sm.dt <- model.data.dt[myx, ]
# 			sm.dt$dataset <- dt
# 	        sm.dt[,1] <- scale(sm.dt[,1])
# 	        sm.dt[,2] <- scale(sm.dt[,2])
# 			sampled.data[[dt]] <- sm.dt
# 			}
# 		sampled.data <- do.call(rbind, sampled.data)

# 	    fixedEff <- tryCatch({
# 	        test.mod <- lmer(y ~ (x + 0| dataset) + x + 0, sampled.data, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
# 	        fixef(test.mod)

# 	        }, error = function(e){ NA_real_})
# 	    t[i] <- fixedEff

# 	}

# 	boot.array <- do.call(cbind,SampleArrayList)

# 	return(list(t = t, boot.array = boot.array))
# }


# ## NOTE, stratified bootstrap modifications not implemented here
# empinf.reg.manualBootArray <- function(t, model.data, boot.array){

# 	R <- length(t)
# 	n <- NROW(model.data)

# 	f <- boot.array

# 	X <- f/n
# 	X <- X[,-1]
# 	beta <- coefficients(glm(t ~ X))[-1L]
# 	l <- numeric(n)

# 	l[-1] <- beta
# 	l <- l - mean(l)
# 	return(l)

# }



# boot.res <- readRDS(file.path(myBootDir, paste0(make.names(drug),"_", make.names(tissue), "_", gene, "_boot_", bootR, "_out.RDS")))

boot.res.unparsed <- readLines(myInFl)


boot.res <- parseBootResults(boot.res.unparsed)

## TODO:: need to write my own code to do regression based 
## imperical influence estimation
## as the current implementation requires me to be able to
## sample based on indicies in dataset
## while I have nested sampling (maybe I can also find a 
## way around this with weights? but unlikely)

# L <- empinf(data=model.data, statistic = )

if(any(is.na(boot.res$t[,1]))){
  warning(paste0("There are: ", sum(is.na(boot.res$t[,1])), " NA results of the bootstrap, equaling ", 100*sum(is.na(boot.res$t[,1]))/length(boot.res$t[,1]), " percent of the data."))
  if(sum(is.na(boot.res$t[,1]))/length(boot.res$t[,1])>0.2){
    stop("Greater than 20 percent of the data is NA.")
  }
}

exp_p <- pnorm(abs(boot.res$t0)/sd(boot.res$t[,1], na.rm=TRUE), lower.tail=FALSE)

optimP <- function(p.val){
  min(abs(boot.ci(boot.res, conf = 1-p.val, type="perc")$perc[4:5]))
}

p.val.boot <- optimise(optimP, c(0,1), tol = min(exp_p/4, 1e-4))$minimum


ci.95.boot <- boot.ci(boot.res, conf = 0.95, type="perc")

boot.sig.res <- list("p.value" = p.val.boot, t0 = boot.res$t0, t1_star = mean(boot.res$t), ci = ci.95.boot, LR = R, bootR = bootR)

saveRDS(boot.sig.res, file=file.path(myPvalDir, make.names(paste("bootSig", gene, drug, tissue, "out.rds", sep="_"))))


