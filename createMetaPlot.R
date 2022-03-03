library(PharmacoGx)
library(meta)

args <- commandArgs(trailingOnly = TRUE)

method <- "perm"

# psetName <- args[1]

drug <- args[1]

tissue <- args[2]

home <- Sys.getenv("HOME")
scratch <- Sys.getenv("SCRATCH")

myDataDir <- file.path(scratch, paste0("pearson_",method,"_res"))
myOutDir <- file.path(scratch, paste0(method, "_meta_plots"))

geneInfo <- read.csv("~/geneInfo.csv")

if(!file.exists(myOutDir)) dir.create(myOutDir)

fls <- list.files(path=file.path(myDataDir), pattern = paste0("*_", drug, "_", tissue, "*"), full.names = TRUE)

sig.list <- sapply(fls, readRDS)

genes_to_check <- unique(unlist(sapply(sig.list, function(x) {
	mygenes <- names(which(1==x[,1,"significant"]))
	mygenes <- gsub(mygenes, pat="\\.[0-9]+$", rep="")
	return(mygenes)
	})))


getTrange <- function(df){
	return(qt(0.975, df) - qt(0.025, df))
}

for(gene in genes_to_check){

	geneCorData <- t(sapply(sig.list, function(x) {
			if(any(grepl(pat=gene, x=rownames(x)))) {
				return(x[grep(pat=gene, x=rownames(x)),,])
			}
			return(x[1,,,drop=FALSE][NA,,])
		}))
	rownames(geneCorData) <- sapply(sig.list, function(x) return(x@PSetName))
	sterr <- (geneCorData[,"upper"] - geneCorData[,"lower"])/getTrange(geneCorData[,"df"])
	# meta.res <- metacor(geneCorData[,"estimate"], n=geneCorData[,"n"], studlab = rownames(geneCorData), hakn = TRUE, adhoc.hakn = "ci", method.tau = "REML")
    meta.res <- metagen(sm="COR", TE = geneCorData[,"estimate"], seTE = sterr, df = geneCorData[,"df"], lower = geneCorData[,"lower"], upper = geneCorData[,"upper"],
    		studlab = rownames(geneCorData), pval = geneCorData[,"pvalue"], approx.TE = "", approx.seTE = "", hakn = TRUE, adhoc.hakn = "ci", method.tau = "REML")
	

    meta.res$pval <- format(meta.res$pval)

    geneName <- geneInfo[grep(gene, geneInfo[,"gene_id"]),"gene_name"]

    pdf(file.path(myOutDir, paste0("nometa_plot_pearson_",method,"_res_", drug, "_", tissue,"_", gene, ".pdf")), width=15, height=5)
	forest(meta.res, comb.fixed=TRUE, comb.random=TRUE, layout="meta", calcwidth.hetstat = TRUE, calcwidth.fixed=TRUE, calcwidth.random=TRUE,
		rightcols=c("effect", "ci", "w.fixed", "w.random", "pval"), text.addline1 = paste("Gene:", unique(geneName)))
	dev.off()

}




### forest(meta.res,  xlim=c(-0.2, 0.7), test.overall.fixed=TRUE, test.overall.random=TRUE,comb.fixed=TRUE, comb.random=TRUE,  scientific.pval=FALSE,  layout="meta", calcwidth.hetstat = TRUE, calcwidth.tests=TRUE, rightcols=c("effect", "ci", "w.fixed", "w.random", "pval"), digits.addcols = 6)



