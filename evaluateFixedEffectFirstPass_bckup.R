
### This file will evaluate the fix effect biomarker results in 2 ways: 
### 1. See whether the known biomarkers are significant
### 2. See if there is enrichment for pathways in biomarkers for a single drug


## Right now, we will be using simple bonferroni correction. This is apt to change to a per-gene correction, or a FDR
library(data.table)

alpha <- 0.05

corP <- 60671

permRes <- fread("biomarker_res/BiomarkerPermRes.csv") 


permRes[, p.adj := PermResP * corP]



permRes[,sum(p.adj < alpha), c("Drug", "Tissue")]


permSig <- permRes[p.adj < alpha]


permSig[,.N, c( "Drug","Tissue")]

gene_info <- fread("geneInfo.csv")

gene_info[,V1 := gsub(V1, pat="\\.[0-9]+$", rep="")]

gene_symbol <- gene_info[,.(V1, gene_name)]

colnames(gene_symbol)[1] <- "Gene"

permSig <- merge(permSig, gene_symbol, by="Gene", all.x=TRUE)





### Lets look at the bootstrap random effects as well here:

bootRes <- fread("biomarker_res/BiomarkerBootRes.csv")

bootRes[, p.adj := BootResP * corP]

bootRes[,sum(p.adj < alpha), c("Drug", "Tissue")]


bootSig <- bootRes[p.adj < alpha]


bootSig <- merge(bootSig, gene_symbol, by="Gene", all.x=TRUE)

allSig <- rbind(permSig[,.(Gene, Drug, Tissue, p.adj, gene_name)], bootSig[,.(Gene, Drug, Tissue, p.adj, gene_name)])



top_markers <- fread("Top biomarkers_Farnoosh.csv", select=1:12, header=TRUE)


top_markers <- top_markers[grep(pat="EXPR", x=Alteration.type),]

top_markers <- top_markers[compound %in% allSig$Drug] 

colnames(top_markers)[1:2] <- c("Drug", "gene_name")

merge(top_markers, allSig, all.x=TRUE, by=c("Drug", "gene_name"))

# That was anti-climatic

## Now lets do a pathway analysis

doFisher <- function(genes, gsc, minSize=1,maxSize=Inf){

	universe <- union(unique(unlist(gsc$gsc)), genes)
	return(fora(gsc$gsc, genes, universe, minSize, maxSize))

}

library(piano)
library(fgsea)

go.bp <- piano::loadGSC("pathways/c5.bp.v7.1.symbols.gmt")

allSig <- allSig[order(Drug, Tissue)]

geneLists <- split(allSig, by=c("Drug","Tissue"))

pathRes <- lapply(geneLists, function(x) return (doFisher(genes=x$gene_name, gsc=go.bp, maxSize=200)))

allPathRes <- rbindlist(pathRes)

allPathRes[,FDR := padj]
allPathRes <- allPathRes[,.(pathway, pathway, pval, padj)]
colnames(allPathRes) <- c("pathway", "Description", "p.Val", "FDR")
fwrite(allPathRes, file="fullGoResCytoScape.txt", sep="\t")

for(name in names(pathRes)){
	drug <- strsplit(name, split="\\.")[[1]][1]
	tissue <- strsplit(name, split="\\.")[[1]][2]
	pathRes[[name]][,Drug := drug]
	pathRes[[name]][,Tissue := tissue]
}
allPathRes <- rbindlist(pathRes)
fwrite(allPathRes, file="fullGoResShiny.txt", sep="\t")


fdrMat <- sapply(pathRes, function(x) return((x[order(pathway),padj])))
rownames(fdrMat) <- pathRes[[1]][order(pathway),pathway]

rowAll0 <- apply(fdrMat, 1, function(x) all(x==1)) 
colAll0 <- apply(fdrMat, 2, function(x) all(x==1)) 

pheatmap(fdrMat[!rowAll0,!colAll0], color = colorRampPalette(rev(brewer.pal(n = 7, name =
       "OrRd")))(100), clustering_method = "ward.D2")

sigAssociations <- apply(fdrMat, 2, function(x) return(x[x<0.1]))

sigAssociations <- sigAssociations[sapply(sigAssociations, length) != 0]

sort(table(unlist(sapply(sigAssociations, names))))

writeLines(names(table(unlist(sapply(sigAssociations, names)))), con="ora_significant_pathways_go_bp.txt")

# ## Using canonical gene sets instead

# canonical <- piano::loadGSC("c2.cp.v7.1.symbols.gmt")

# allSig <- allSig[order(Drug, Tissue)]

# geneLists <- split(allSig, by=c("Drug","Tissue"))

# pathRes <- lapply(geneLists, function(x) return (doFisher(genes=x$gene_name, gsc=canonical)))

# fdrMat <- sapply(pathRes, function(x) return((x[order(pathway),padj])))
# rownames(fdrMat) <- pathRes[[1]][order(pathway),pathway]

# rowAll0 <- apply(fdrMat, 1, function(x) all(x==1)) 
# colAll0 <- apply(fdrMat, 2, function(x) all(x==1)) 

# pheatmap(fdrMat[!rowAll0,!colAll0], color = colorRampPalette(rev(brewer.pal(n = 7, name =
#        "OrRd")))(100), clustering_method = "ward.D2")

# sigAssociations <- apply(fdrMat, 2, function(x) return(x[x<0.1]))

# sigAssociations[sapply(sigAssociations, length) != 0]


# writeLines(names(table(unlist(sapply(sigAssociations, names)))), con="ora_significant_pathways_c2.txt")


### Out of curiousity, we will do a network analysis to calculate mean distance between markers and targets

#### start with GeneMania, because I like the integrative network

library(Matrix)
library(igraph)

geneManiaNetwork <- fread("GeneMania/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt")

geneLevels <- geneManiaNetwork[,union(unique(Gene_A), unique(Gene_B))]
geneManiaNetwork[,Gene_A := factor(Gene_A, levels = geneLevels)]
geneManiaNetwork[,Gene_B := factor(Gene_B, levels = geneLevels)]


geneManiaSM <- sparseMatrix(i=as.numeric(geneManiaNetwork$Gene_A), j = as.numeric(geneManiaNetwork$Gene_B), x=geneManiaNetwork[,Weight], symmetric=TRUE)

geneManiaGraph <- graph_from_adjacency_matrix(geneManiaSM, weighted=TRUE, mode="undirected")

rownames(geneManiaSM) <- colnames(geneManiaSM) <- geneLevels

drugTargetInfo <- fread("~/Code/Github/pachyderm/Annotations/DrugTargetCompilation_updated.csv")


drugTargetInfo <- drugTargetInfo[BHKLAB.DRUGID %in% allSig$Drug]

drugTargetInfo <- drugTargetInfo[,.("TARGET"= unique(TARGET_NAME)),BHKLAB.DRUGID]

drugTargetInfo <- drugTargetInfo[complete.cases(drugTargetInfo),]

drugTargetInfo <- drugTargetInfo[drugTargetInfo$TARGET %in% gene_info$gene_name,]

drugTargetInfo[,TARGET := gene_info$V1[match(TARGET, gene_info$gene_name)]]

distance.list <- list()

for(drug in drugTargetInfo[,unique(BHKLAB.DRUGID)]){

	targets <- drugTargetInfo[BHKLAB.DRUGID == drug, TARGET]

	markers <- allSig[Drug==drug, Gene]

	tg <- unique(na.omit(match(targets, colnames(geneManiaSM))))

	mk <- unique(na.omit(match(markers, colnames(geneManiaSM))))


	distance.list[[drug]] <- distances(geneManiaGraph, tg, mk)

}

tg_random <- sample(length(geneLevels), 30)
mk_random <- sample(length(geneLevels), 30)

random_dists <- distances(geneManiaGraph, tg_random, mk_random)

# graph_distances <- distance_table(geneManiaGraph) # probably way too long to run


random_dists.m <- reshape2::melt(random_dists)
distance.list.m <- reshape2::melt(unlist(distance.list))


library(ggplot2)

toPlot <- data.frame(Cat = c(rep("Random" , times=length(random_dists.m[,3])), rep("Marker to Target", times=length(distance.list.m[,1]))), "Distance" = c(random_dists.m[,3], distance.list.m[,1]))


ggplot(toPlot, aes(x=log10(Distance), fill=Cat)) + geom_density(alpha=0.6) + scale_fill_manual(	values=c("red", "royalblue"))

### Should try this with a simple ORA analysis. 




