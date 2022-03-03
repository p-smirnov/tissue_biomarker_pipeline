## merged to runFixedEffectPerm


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

library(RhpcBLASctl)
print(RhpcBLASctl::blas_get_num_procs())
# RhpcBLASctl::blas_set_num_threads(1)
# registerDoParallel(40)

bootR <- "1e.07"

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
R <- as.numeric(args[4]) ## In this case, its the number of samples to use for bootstrap

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")
project <- Sys.getenv("PROJECT")


myInFl <- args[5]
myPvalDir <- file.path(project, paste0(method, "_meta_perm_sig"))


if(!file.exists(myPvalDir)) dir.create(myPvalDir)



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


bootCor <- function(model.data,R){

    N <- nrow(model.data)

    t <- unlist(mclapply(seq_len(R), function(i, model.data) {

      model.data.sampled <- getBootSample(model.data)
      denom  <- sum(model.data.sampled[,"x"]^2)
      xt <- model.data.sampled[,"x"]/denom


      crossprod(xt,model.data.sampled[,"y"])[1]
    }, model.data = model.data), mc.cores=nthread)


     model.data.scaled <- do.call(rbind, lapply(model.data$dataset, function(dt){
        sm.dt <- model.data[model.data$dataset == dt,]
        sm.dt[,1] <- scale(sm.dt[,1])
        sm.dt[,2] <- scale(sm.dt[,2])
        return(sm.dt)
        }))

    denom  <- sum(model.data.scaled[,"x"]^2)
    xt <- model.data.scaled[,"x"]/denom

    t0 <- crossprod(xt, model.data.scaled[,"y"])[1]


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



### First we calculate the p-value and do a bootstrap for the 95% confidence interval estimation

perm.res <- readRDS(myInFl)

p.value <- mean(abs(perm.res$t0) <= abs(perm.res$t))
if(p.value == 0){
  p.value <- 1/(length(perm.res$t)+1)
}

model.data <- data.frame(x=perm.res$data.x, y=perm.res$data.y,dataset=perm.res$data.dataset)

boot.res <- bootCor(model.data, R) 


ci.95.boot <- boot.ci(boot.res, conf = 0.95, type="perc")

perm.sig.res <- list("p.value" = p.value, ci = ci.95.boot, LR = R, R = length(perm.res$t))

saveRDS(perm.sig.res, file=file.path(myPvalDir, make.names(paste("permSig", gene, drug, tissue, "out.rds", sep="_"))))

