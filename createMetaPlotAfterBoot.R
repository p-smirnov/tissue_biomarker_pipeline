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

### Now lets plot our results
R <- as.numeric(args[4]) ## In this case, its the number of samples to use for L estimation

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- Sys.getenv("DATA")
mySigDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myBootDir <- file.path(scratch, paste0(method, "_meta_bootstrap"))
myFigDir <- file.path(scratch, paste0(method, "_meta_boot_figures"))
myPvalDir <- file.path(scratch, paste0(method, "_meta_boot_sig"))


if(!file.exists(myFigDir)) dir.create(myFigDir)
if(!file.exists(myPvalDir)) dir.create(myPvalDir)

myBootDir <- file.path(scratch, "juliaBoot")

bootR <- "10000000"

boot.sig.res <- readRDS(file=file.path(myPvalDir, make.names(paste("bootSig",drug, tissue, gene, "out.rds", sep="_"))))



toRunExtended <- read.csv("~/tissueDrugToRunExtended.txt", header=FALSE)
pSets <- toRunExtended[toRunExtended[,2] == drug & toRunExtended[,3] == tissue, 1]

geneInfo <- read.csv("~/geneInfo.csv")

fls <- paste0(file.path(mySigDir,paste("signature", make.names(pSets), make.names(drug), make.names(tissue), sep="_")), ".rds")
sig.list <- sapply(fls, readRDS)


getTrange <- function(df){
    return(qt(0.975, df) - qt(0.025, df))
}


geneCorData <- t(sapply(sig.list, function(x) {
        if(any(grepl(pat=gene, x=rownames(x)))) {
            return(x[grep(pat=gene, x=rownames(x)),,])
        }
        return(x[1,,,drop=FALSE][NA,,])
    }))



rownames(geneCorData) <- sapply(sig.list, function(x) return(x@PSetName))
sterr <- (geneCorData[,"upper"] - geneCorData[,"lower"])/getTrange(geneCorData[,"df"])


library(ggplot2)




### Adapted from: 
### https://gist.github.com/webbedfeet/7031404fc3f500f6258e
credplot.gg <- function(d, xlab='Variable', ylab="Y", main="Plot", meta.label="Meta"){
  # d is a data frame 
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  # d$annot labels the graph
  if(is.character(d$y)) d$y <- factor(d$y, levels=rev(d$y))
  if(meta.label %in% d$y){
    d.meta <- d[d$y == meta.label,]
    d.size <- d.meta[,"size"]
    d.meta <- rbind(d.meta, d.meta, d.meta, d.meta)
    d.meta$y <- as.numeric(d.meta$y)
    d.meta[4,c("x","y")] <-  as.numeric(d.meta[1,c("xlo","y")])
    d.meta[3,c("x","y")] <-  as.numeric(d.meta[1,c("x","y")]) + c(0,0.1)*d.size
    d.meta[2,c("x","y")] <-  as.numeric(d.meta[1,c("xhi","y")])
    d.meta[1,c("x","y")] <-  as.numeric(d.meta[1,c("x","y")]) - c(0,0.1)*d.size
    d[d$y == meta.label,c("xlo", "xhi")] <- NA_real_
    # d.meta$xhi <- NA_real_
  }
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, xmin=xlo, xmax=xhi, label=annot))+
    geom_errorbarh(height = 0.15, show.legend=FALSE)+ geom_point(size=d$size) + 
    geom_vline(xintercept = 0, linetype=2,show.legend=FALSE) + xlim(c(-1,1)) + 
    coord_cartesian(clip = 'off') +
    xlab(xlab) +
    ylab(ylab) + theme_bw() + geom_text(aes(x=Inf, y=y), nudge_x = 1.0, size=4, hjust="outward",show.legend=FALSE)+# scale_y_discrete(breaks=dd$y)  + 
    theme(plot.margin = unit(c(1,8,1,1), "lines"),legend.position="bottom") + ggtitle(main) 
  
  if(meta.label %in% d$y){
    p <- p + geom_polygon(data=d.meta)
  }
  return(p)
}


geneName <- geneInfo[grep(gene, geneInfo[,"gene_id"]),"gene_name"]


dd <- data.frame(y=rownames(geneCorData), x=geneCorData[,"estimate"],
                 xlo = geneCorData[,"lower"], xhi = geneCorData[,"upper"], 
                 annot = geneCorData[,"pvalue"], size = log2(geneCorData[,"n"]))
dd <- rbind(dd, data.frame(y = "Meta", x = mean(boot.res$t[,1]), xlo = boot.sig.res$ci$bca[4],
                  xhi = boot.sig.res$ci$bca[5], annot = boot.sig.res$p.value, size = 2))
dd$annot <- paste0("  p value: ", format(dd$annot, scientific=T, digits=3))


pdf(file.path(myFigDir, paste0("plot_pearson_",method,"_res_", make.names(drug), "_", make.names(tissue),"_", gene, "_pval_", format(boot.sig.res$p.value, scientific=T, digits=3),".pdf")), width=15, height=5)
# gtbl <- tableGrob(data.frame(p.value = dd$annot, N = c(geneCorData[,"n"], NA_real_)), theme=theme_minimal())
# gtbl$heights <- unit(rep(1/nrow(gtbl), nrow(gtbl)), "npc")
p <- credplot.gg(dd, "Estimate","Dataset", main = paste0("Drug: ", drug, " Gene: ", geneName)) ## Change to gene name
print(p)
# grid.arrange(p, gtbl, ncol=2)
dev.off()



### forest(meta.res,  xlim=c(-0.2, 0.7), test.overall.fixed=TRUE, test.overall.random=TRUE,comb.fixed=TRUE, comb.random=TRUE,  scientific.pval=FALSE,  layout="meta", calcwidth.hetstat = TRUE, calcwidth.tests=TRUE, rightcols=c("effect", "ci", "w.fixed", "w.random", "pval"), digits.addcols = 6)



