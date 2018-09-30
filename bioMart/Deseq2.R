library(DESeq2)
library(limma)
library(pasilla)
library(SRAdb)

data_geo <- readRDS('./data/data_geo.rds')
data_geo <- read.table('./data/GSE64912_Read_counts 2.txt',
                       header = T, row.names = 1)
data("pasillaGenes")
# 表达矩阵
exprSet <- counts(pasillaGenes)
exprSet2 <- as.matrix(data_geo)

# 分组因子
group_list <- pasillaGenes$condition
pheno <- phenoData(data_geo)
group_list2 <- unname(sapply(as.character(pheno@data$characteristics_ch1),
                      function(x) unlist(strsplit(x, split = ': '))[2]))

# 分组矩阵？
col_data <- data.frame(row.names = colnames(exprSet), group_list = group_list)
col_data2 <- data.frame(group_list = group_list2, row.names = colnames(exprSet2))

dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = col_data,
                              design = ~ group_list)
dds2 <- DESeqDataSetFromMatrix(countData = exprSet2, colData = col_data2,
                               design = ~ group_list2)

dds_2 <- DESeq(dds)
resultsNames(dds_2)
res <- results(dds_2, contrast = c('group_list', 'treated', 'untreated'))
resOrdered <- res[order(res$padj),]
resOrdered <- as.data.frame(resOrdered)


# SRA ---------------------------------------------------------------------

sqlfile <-'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
sra_con <- dbConnect(SQLite(), sqlfile)