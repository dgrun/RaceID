## install required packages (only at first time)
##install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit"))

## load class definition and functions
source("RaceID_class.R")

## input data
x <- read.csv("transcript_counts_intestine.xls",sep="\t",header=TRUE)
rownames(x) <- x$GENEID
# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]

## RaceID
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
# k-means clustering
sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000)
# compute t-SNE map
sc <- comptsne(sc,rseed=15555)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)

## diagnostic plots
# gap statistics
plotgap(sc)
# silhouette of k-means clusters
plotsilhouette(sc)
# Jaccard's similarity of k-means clusters
plotjaccard(sc)
# barchart of outlier probabilities
plotoutlierprobs(sc)
# regression of background model
plotbackground(sc)
# dependence of outlier number on probability threshold (probthr)
plotsensitivity(sc)
# heatmap of k-means cluster
clustheatmap(sc,final=FALSE,hmethod="single")
# heatmap of final cluster
clustheatmap(sc,final=TRUE,hmethod="single")
# highlight k-means clusters in t-SNE map
plottsne(sc,final=FALSE)
# highlight final clusters in t-SNE map
plottsne(sc,final=TRUE)
# highlight cell labels in t-SNE map
plotlabelstsne(sc,labels=sub("(\\_\\d+)","",names(sc@ndata)))
# highlight groups of cells by symbols in t-SNE map
plotsymbolstsne(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))
# highlight transcirpt counts of a set of genes in t-SNE map, e. g. all Apoa genes
g <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", "Apoa5__chr9")
plotexptsne(sc,g,n="Apoa genes",logsc=TRUE)

## identification of marker genes
# differentially regulated genes in each cluster compared to the full ensemble
cdiff <- clustdiffgenes(sc,pvalue=.01)

## write results to text files
# final clusters 
x <- data.frame(CELLID=names(sc@cpart),cluster=sc@cpart)
write.table(x[order(x$cluster,decreasing=FALSE),],"cell_clust.xls",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
# differentially expressed genes in cluster
for ( n in names(cdiff) ) write.table(data.frame(GENEID=rownames(cdiff[[n]]),cdiff[[n]]),paste(paste("cell_clust_diff_genes",sub("\\.","\\_",n),sep="_"),".xls",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

# differentially expressed genes between two sets of clusters, e. g. cluster 1 and clusters 2,3
d <- diffgenes(sc,cl1=1,cl2=c(2,3),mincount=5)
plotdiffgenes(d,gene=names(d$z)[1])


