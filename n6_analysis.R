########## Installing packages #################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

install.packages("tidyverse")
install.packages("dplyr")
install.packages("plyr")
install.packages("tidyr")
install.packages("factoextra")
install.packages("pheatmap")
########### Loading packages #################
library(plyr)
library(tidyr)
library(dplyr)
library(readxl)
library(DESeq2)
library("org.Hs.eg.db")
library(ggplot2)
library(ComplexHeatmap)
library(factoextra)
library("pheatmap")
###############################################
#set working directory in the format /path/to/working/directory
setwd("C:/Users/chidi/Documents/Summer_2023/")

data_n6 <- read_excel("C:/Users/chidi/Documents/Summer_2023/n6_tissue_rawread.xls") 

#converting to dataframe
data_n6 <- data.frame(data_n6)

#setting row names
rownames(data_n6) <- data_n6$EnsemblID

#deleting GeneID and EnsembleID columns 
data_n6 <- data_n6[,-1]

#input correct colnames
colnames(data_n6) <-c("11B_A", "11B_C", "11B_D", "11B_ML",	"11B_MU",
"11B_P",	"13B_A",	"13B_C",	"13B_D",	"13B_ML",	"13B_MU",
"13B_P",	"29B_A", "29B_C", "29B_D", "29B_ML", "29B_MU",	"29B_P", "31B_A", "31B_C", "31B_D",
"31B_ML", "31B_MU", "31B_P", "34B_A",	"34B_C",	"34B_D",	"34B_ML", "34B_MU", "34B_P",
"35B_A",	"35B_C", "35B_D",	"35B_ML", "35B_MU", "35B_P",	"37B_A",	"37B_C",	"37B_D",	
"37B_ML",	"37B_MU",	"37B_P",	"40B_A",	"40B_C",	"40B_D",	"40B_ML",	"40B_MU",	
"40B_P",	"44B_A",	"44B_C",	"44B_D",	"44B_ML",	"44B_MU",	"44B_P",	"45B_A",	"45B_C",	
"45B_D",	"45B_ML",	"45B_MU",	"45B_P",	"46B_A",	"46B_C",	"46B_D",	"46B_ML",	"46B_MU",
"46B_P",	"52A_A",	"52A_C",	"52A_D",	"52A_ML",	"52A_MU",	"52A_P",	"53A_A",	"53A_C",	
"53A_D",	"53A_ML",	"53A_MU",	"53A_P",	"55A_A",	"55A_C",	"55A_D",	"55A_ML",	"55A_MU",	
"55A_P",	"56A_A",	"56A_C",	"56A_ML",	"56A_MU",	"56A_P",	"56_D",	"57A_A",	"57A_C",	"57A_D",	
"57A_ML",	"57A_MU",	"57A_P",	"59A_A",	"59A_C",	"59A_D",	"59A_ML",	"59A_MU",	
"59A_P",	"61A_A",	"61A_C",	"61A_D",	"61A_ML",	"61A_MU",	"61A_P",	"62A_A",	"62A_C",	
"62A_D",	"62A_ML",	"62A_MU",	"62A_P",	"9B_A", "9B_C",	"9B_D",	"9B_ML",	"9B_MU",	"9B_P"

)

#reading in metadata
metaData_n6 = read.delim("metaData_n6.txt", header=TRUE)

#add sample name to rowname and delete sample column
rownames(metaData_n6) <- metaData_n6$Sample

#check to see whether row names in metaData match column names in dataFrame
all(colnames(data_n6) %in% rownames(metaData_n6))

#the values are also in the same order? 
all(colnames(data_n6) == rownames(metaData_n6))

#constructing a DeSeq2 object
dds_n6 <-DESeqDataSetFromMatrix(countData=data_n6, 
                                 colData=metaData_n6,
                                 design= ~Location)

#specify factor level
dds_n6$Location <- relevel(dds_n6$Location, ref="Myo_Lower")

levels(dds_n6$Location)

#checking model matrix
model.matrix(~Location, data = metaData_n6)

#pre-filtering: keeping only rows with at least 10 reads
keep <- rowSums(counts(dds_n6)) >= 10
dds_n6 <- dds_n6[keep,]

dds_n6$group <- factor(paste0(dds_n6$Location, dds_n6$Status))
design(dds_n6) <- ~ group
dds_n6$group <- relevel(dds_n6$group, ref="Myo_LowerTL")

levels(dds_n6$group)

dds_n6 <- DESeq(dds_n6)

levels(dds_n6$group)

dds_n6 <- DESeq(dds_n6)
resultsNames(dds_n6)
results(dds_n6, contrast=c("group", "IB", "IA"))

res_n6_AvML<- results(dds_n6, contrast=c("group","", ""))

res_n6_AvML<- results(dds_n6, contrast=c("Location","Amnion", "Myo_Lower"))
res_n6_CvML<- results(dds_n6, contrast=c("Location","Chorion", "Myo_Lower"))
res_n6_DvML<- results(dds_n6, contrast=c("Location","Decidua", "Myo_Lower"))
res_n6_MUvML<-results(dds_n6, contrast=c("Location","Myo_Upper", "Myo_Lower"))
res_n6_PvML<- results(dds_n6, contrast=c("Location","Placenta", "Myo_Lower"))

#sort results by p-value
res_n6_AvML <- res_n6_AvML[order(res_n6_AvML$padj),]
head(res_n6_AvML)

res_n6_CvML <- res_n6_CvML[order(res_n6_CvML$padj),]
head(res_n6_CvML)

res_n6_DvML <- res_n6_DvML[order(res_n6_DvML$padj),]
head(res_n6_DvML)

res_n6_MUvML <- res_n6_MUvML[order(res_n6_MUvML$padj),]
head(res_n6_MUvML)

res_n6_PvML <- res_n6_PvML[order(res_n6_PvML$padj),]
head(res_n6_PvML)

#convert result of DESeq2 to a dataframe
resDF_n6_AvML <- as.data.frame(res_n6_AvML)
resDF_n6_CvML <- as.data.frame(res_n6_CvML)
resDF_n6_DvML <- as.data.frame(res_n6_DvML)
resDF_n6_MUvML <- as.data.frame(res_n6_MUvML)
resDF_n6_PvML <- as.data.frame(res_n6_PvML)

#adding gene IDs
resDF_n6_AvML$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(resDF_n6_AvML), 
                              keytype = "ENSEMBL", column = "SYMBOL")

resDF_n6_CvML$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(resDF_n6_CvML), 
                               keytype = "ENSEMBL", column = "SYMBOL")

resDF_n6_DvML$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(resDF_n6_DvML), 
                               keytype = "ENSEMBL", column = "SYMBOL")

resDF_n6_MUvML$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(resDF_n6_MUvML), 
                               keytype = "ENSEMBL", column = "SYMBOL")

resDF_n6_PvML$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(resDF_n6_PvML), 
                               keytype = "ENSEMBL", column = "SYMBOL")

#removing NAs value in padj, filtering for DEG, printing to cvs
DEG_n6_AvML <- resDF_n6_AvML[(abs(resDF_n6_AvML$log2FoldChange) > 1.5) 
                           & (resDF_n6_AvML$padj) < 0.05 , ]

DEG_n6_CvML <- resDF_n6_CvML[(abs(resDF_n6_CvML$log2FoldChange) > 1.5) 
                             & (resDF_n6_CvML$padj) < 0.05 , ]

DEG_n6_DvML <- resDF_n6_DvML[(abs(resDF_n6_DvML$log2FoldChange) > 1.5) 
                             & (resDF_n6_DvML$padj) < 0.05 , ]

DEG_n6_MUvML <- resDF_n6_MUvML[(abs(resDF_n6_MUvML$log2FoldChange) > 1.5) 
                             & (resDF_n6_MUvML$padj) < 0.05 , ]

DEG_n6_PvML <- resDF_n6_PvML[(abs(resDF_n6_PvML$log2FoldChange) > 1.5) 
                               & (resDF_n6_PvML$padj) < 0.05 , ]

write.csv(DEG_n6_AvML, file = "topGenes_n6_AvML.csv")
write.csv(DEG_n6_CvML, file = "topGenes_n6_CvML.csv")
write.csv(DEG_n6_DvML, file = "topGenes_n6_DvML.csv")
write.csv(DEG_n6_MUvML, file = "topGenes_n6_MUvML.csv")
write.csv(DEG_n6_PvML, file = "topGenes_n6_PvML.csv")

#RHOQ 
RHOQ <- plotCounts(dds_n6, "ENSG00000119729", intgroup = c("group"), returnData = TRUE) 
ggplot(RHOQ, aes(x=group, y = count)) + 
  geom_point(aes(x=group), colour = "green", size = 3) + 
  ggtitle("RHOQ expression") +
  theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1))+
  xlab("Group") + ylab("normalized count")

#RGS2
plotCounts(dds_n6, "ENSG00000116741", intgroup = c("Location")) 

#PDE4D
plotCounts(dds_n6, "ENSG00000113448", intgroup = c("Location")) 

#MYOCD
plotCounts(dds_n6, "ENSG00000141052", intgroup = c("Location")) 
"ENSG00000141052" %in% rownames(dds_n6)

#EP1
plotCounts(dds_n6, "ENSG00000125384", intgroup = c("Location")) 
"ENSG00000125384" %in% rownames(dds_n6)

#EP2
plotCounts(dds_n6, "ENSG00000113448", intgroup = c("Location")) 

#EP3
plotCounts(dds_n6, "ENSG00000050628", intgroup = c("Location")) 

#EP4
plotCounts(dds_n6, "ENSG00000171522", intgroup = c("Location")) 

################ PCA ###############################
rld_n6 <- vst(dds_n6)
head(assay(rld_n6))

#PCA plot using log transformed data

plotPCA(rld_n6, intgroup = "group")

rld_mat_n6 <- assay(rld_n6)
pca_n6 <- prcomp(t(rld_mat_n6))
df_n6 <- cbind(metaData_n6, pca_n6$x)
ggplot(df_n6) + geom_point(aes(x=PC1, y=PC2)) + 
  geom_text(aes(x=PC1, y=PC2, label=rownames(df_n6)), vjust = 2)

############### Heatmap ##########################

select <- order(rowMeans(counts(dds_n6,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_n6)[,c("Location","Status")])

df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(select), 
                       keytype = "ENSEMBL", column = "SYMBOL")
pheatmap(df)

pheatmap(assay(dds_n6)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)



