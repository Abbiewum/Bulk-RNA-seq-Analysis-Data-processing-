getwd()
setwd("C:\\Users\\abbie\\Desktop\\counts")

if(!require("BiocManager",quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install(version="3.19")

BiocManager::available()

BiocManager::install("DESeq2")
library("DESeq2")


install.packages("rmarkdown")
library("rmarkdown")

counts_SRR28420795 <- read.table("counts_SRR28420795.txt", header=TRUE, row.names=1)
counts_SRR28420796 <- read.table("counts_SRR28420796.txt", header=TRUE, row.names=1)
counts_SRR28420797 <- read.table("counts_SRR28420797.txt", header=TRUE, row.names=1)
counts_SRR28420798 <- read.table("counts_SRR28420798.txt", header=TRUE, row.names=1)

View(counts_SRR28420795)


counted_data <- cbind(counts_SRR28420795, counts_SRR28420796, counts_SRR28420797, counts_SRR28420798)


counted_data <- counted_data[, seq(6, ncol(counted_data), by = 6)]


colnames(counted_data) <- gsub("^X\\.home\\.abby\\.|\\.sorted\\.bam$", "", colnames(counted_data))


View(counted_data)

counted_data <- as.matrix(counted_data)

head(counted_data)

condition = c("Treated", "Control", "Treated", "Control")

mycols <- c("Control" = "red", "Treated" = "blue")

print(condition)

library(DESeq2)

(coldata <- data.frame(row.names=colnames(counted_data), condition))

dds <- DESeqDataSetFromMatrix(countData=counted_data, colData=coldata, design=~condition)

dds

dds <- DESeq(dds)

png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()


rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

sampleDists <- as.matrix(dist(t(assay(rld))))
install.packages("gplots")
library(gplots)

png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "green", "yellow"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(15, 15), main="Sample Distance Matrix")
dev.off()

# Principal Component Analysis (PCA)

png("PCA.png", w=1000, h=1000)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()


library(ggplot2)

BiocManager::install("genefilter")
library(genefilter)

BiocManager::install("grDevices")
library(grDevices)

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(dds)[select,]))
scores <- data.frame(pc$x, condition)

View(scores)

png("PCA1.png", w=1000, h=1000)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1") 
  + geom_text(aes(label = rownames(scores)), size = 4, vjust = -1.5)
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))
dev.off()


select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
png("Heatmap.png", w=1000, h=1000, pointsize=20)
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(dds)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.off()

res <- results(dds)
table(res$padj<0.05)

res <- res[order(res$padj), ]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

write.csv(resdata, file="diffexpr-results.csv")

hist(res$pvalue, breaks=50, col="grey")


install.packages("calibrate")
library(calibrate)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()



volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
