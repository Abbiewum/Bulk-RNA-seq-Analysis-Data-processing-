BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("biomartr")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")


install.packages("fastmap")
packageVersion("fastmap")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HDO.db")


library("biomaRt")
library("clusterProfiler")
library("enrichplot")
library("biomartr")
library("tidyverse")
library("DOSE")

BiocManager::install("org.Hs.eg.db")

suppressPackageStartupMessages(library("org.Hs.eg.db"))

resdata <- read.csv("diffexpr-results.csv")


genes <- biomartr::organismBM(organism = "Homo sapiens")
genes
View(genes) 


hsapiens_attributes = 
  biomartr::organismAttributes("Homo sapiens") %>% 
  filter(dataset == "hsapiens_gene_ensembl")
hsapiens_attributes
View(hsapiens_attributes)



filters <- biomartr::getFilters(mart    = "ENSEMBL_MART_ENSEMBL", 
                                dataset = "hsapiens_gene_ensembl")
View(filters)


gene_set <- resdata$Gene
result_BM <- biomartr::biomart(
  genes = gene_set,
  mart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id"
)
head(result_BM)



result_df <- data.frame(result_BM)


colnames(result_df) <- c("Gene", "HGNC_symbol", "Entrez_gene_id")


geneID <- result_df$Entrez_gene_id


ego <- enrichGO(
  gene = geneID,
  OrgDb = org.Hs.eg.db,  # Human genome annotation database
  keyType = "ENTREZID",  # Use Entrez gene IDs as keys
  ont = "BP",            # Biological Process ontology. Can either be BP, CC or MF
  pAdjustMethod = "BH",  #  
  pvalueCutoff = 0.05,   # P-value cutoff for significance
  qvalueCutoff = 0.05,   # Adjusted p-value (FDR) cutoff for significance
  readable= TRUE)
ego


png("GO_BP.png", w=1000, h=1000, pointsize=20)
dotplot(ego, showCategory = 15) + ggtitle("Gene Enrichment Plot") + theme_bw() # Show top 15 enriched terms in a dotplot
dev.off()


kegg_enrich <- enrichKEGG(
  gene = geneID,
  organism = "hsa",       # Specify organism code for Homo sapiens
  pvalueCutoff = 0.05,    # P-value cutoff for significance
  qvalueCutoff = 0.05,     # Adjusted p-value (FDR) cutoff for significance
)
View(kegg_enrich)

str(kegg_enrich)
head(kegg_enrich)

install.packages("clusterProfiler")
library("clusterProfiler")

data(geneList, package = "DOSE")
de <- names(geneList)[abs(geneList) > 2]

kegg_enrich <- enrichKEGG(gene = de, organism = 'hsa')

if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
  stop("No enriched pathways found. Check your enrichment results.")
}

library(clusterProfiler)


png("KEEG.png", w=1000, h=1000, pointsize=20)
dotplot(kegg_enrich, showCategory = 10, color = "qvalue", size = "Count")  # Show top 10 enriched pathways in a dotplot
dev.off()






