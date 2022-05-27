## This function makes sure that AACs are recalculated over a common concentration range across cell lines, thereby making them more comparable.
## It does this by trying to maximize the number of measured points remaining, removing points that fall outside the chosen concentration 
## range if measurements for a curve exist past that range, and removing full curves if the whole range was not measured.
## In practice, I find it has been removing about 5% of the data, but that is of course dependent of the dataset. 
standardizeRawDataConcRange <- function(sens.info, sens.raw){
  unq.drugs <- unique(sens.info$drugid)
  
  conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
  conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  conc.ranges[,Var1 := NULL]
  conc.ranges <- conc.ranges[,unique(.SD), drugid]	
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
  l = sq[seq(1,length(sq)-1)];
  r = sq[seq(2,length(sq))];
  .(l=l,r=r)}, drugid]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N>0)
    unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl){
    if(length(indicies) == 1) return(TRUE)
    sq <- seq(indicies[1], indicies[length(indicies)]-1)
    all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][, N], returnConsInts)
  
  names(per.drug.range.indicies) <- conc.ranges.disj[, unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
  
  
  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug) {
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
    
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug) {
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify=FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
  
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)

  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug){
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
      conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })

  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
    return(exp.out.range)
  }))
  
  sens.raw[removed.experiments,,] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]
  
  for(drug in unq.drugs){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    myx <- conc.ranges.kept[drugid==drug,Var1]
    doses <- sens.raw[myx, ,"Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    
    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE
  
  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}


#filter noisy curves from PSet (modified from PharmacoGx, to take into account standardized conc range)
filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
  acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
    #for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
    if (!all(is.na(drug.responses))){
      
      
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
      
      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
      
      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
      
      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
          (delta.sum < epsilon) &
          (max.cum.sum < (2 * epsilon)) &
          (mean(drug.responses$Viability) < mean.viablity)) {
        return (xp)
      }
    }
    
  }, mc.cores=nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  if (is.null(nrow(tt))){
    tt <- matrix(tt, ncol = 2)
  }
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}   

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
        if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}
  
  