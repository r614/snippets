# open a new R script
rm(list = ls())


# load required libraries -------------------------------------------------
#  Read sample gene counts with a csv file.
suppressMessages(library(tidyverse))
suppressMessages(library( "DESeq2" ))
suppressMessages(library("pheatmap"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library(pathview))
suppressMessages(library(gage))
suppressMessages(library(gageData))


# Load data ---------------------------------------------------------------
#  code here....

countdata <- read.csv("countdata.csv")
metadata <- read.csv("metaData.csv")
#Set Row Names
countdata <- countdata %>% remove_rownames %>% column_to_rownames(var="X")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="X")

# Q1 ----------------------------------------------------------------------
#  Write the code for Q1

# Checks number of rows with all zero values.
 
#logdata <- log(select(countdata))
#logdata$avg <- rowMeans(logdata)
#dropped <- nrow(logdata[is.infinite(rowSums(logdata)),])
#logdata <- logdata[is.finite(rowSums(logdata)),]

#cols <- colnames(logdata)
#print(cols)
#cols <- cols[cols != "avg"]

#cat("Number of genes with zero counts in one or more samples:", dropped)

#scaledcounts <- countdata

#for(col in cols){
#  logdata[[col]] <- logdata[[col]] - logdata[["avg"]]
#  scaleVal <- exp(median(logdata[[col]]))
#  print(scaleVal)
#  scaledcounts[[col]] <- scaleVal*scaledcounts[[col]]
#}

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata, design = ~ strain)
dds <- DESeq(dds)
res <- results(dds)


# Q2 ----------------------------------------------------------------------
#  Write the code for Q2

sum(res$padj < 0.1, na.rm=TRUE)



# Q3 --------------------------------------------------------------------

res07 <- results(dds, alpha=0.07)

top_up   <- res07 %>% as.data.frame() %>% na.omit() %>% filter(log2FoldChange > 0, padj < 0.07)
top_down <- res07 %>% as.data.frame() %>% na.omit() %>%filter(log2FoldChange < 0, padj < 0.07) 

total_FDR_cutoff = nrow(top_up) + nrow(top_down)

cat("Number of genes up-regulated:", nrow(top_up), "\n",
    "Number of genes down-regulated:",  nrow(top_down), "\n",
    "Total number of genes at FDR < 0.07:", total_FDR_cutoff, "\n")

# Q4 --------------------------------------------------------------------

top_up_LFC  <- top_up[order(top_up$log2FoldChange), ]
top_up_LFC  <- top_up_LFC[1:20, ]
top_down_LFC <- top_down[order(top_down$log2FoldChange), ]
top_down_LFC<- top_down_LFC[1:20, ]

top_up_padj <- top_up_LFC[order(top_up_LFC$padj), ]
top_down_padj <- top_down_LFC[order(top_down_LFC$padj), ]


cat("Top 20 up-regulated genes:\n\n"); top_up_padj
cat("Top 20 down-regulated genes:\n\n"); top_down_padj



# Q5 --------------------------------------------------------------------

plotMA(res, ylim=c(-2,2))

# Q6 --------------------------------------------------------------------


plotCounts(dds, gene=which.min(res$pvalue), intgroup = "strain")

meanp <- mean(res$pvalue, na.rm=TRUE)

plotCounts(dds, gene=which.min(abs(res$pvalue-meanp)), intgroup = "strain")


# Q7 ----------------------------------------------------------------------

maxIndex <- which.max(res$pvalue)
gene <- rownames(res[maxIndex,])[1]
cat("Gene with max pval is", gene)

high_pval_count <- counts(dds[gene,],)


# Q8 ----------------------------------------------------------------------

res$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res$name <- mapIds(org.Mm.eg.db, keys=row.names(res), column="GENENAME", keytype="ENSEMBL", multiVals="first")

write.csv(res,'resultPath.csv')

# Q9 ----------------------------------------------------------------------


data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
''
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez

gageres <- gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)

pathways <- data.frame(id=rownames(gageres$less), gageres$less) %>% tbl_df() %>% filter(row_number()<=6) %>% .$id %>%  as.character()

ids <- substring(pathways, first = 1, last = 9)

cat("6 downregulated genes are:", ids)



# Q10 ---------------------------------------------------------------------

ntd <- vst(dds)
select <- order(gageres$greater,
                decreasing=TRUE)[1:15]
genes <- order(res$log2FoldChange, decreasing=TRUE)[1:15]

df <- as.data.frame(colData(dds)['strain'])

pheatmap(assay(ntd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
