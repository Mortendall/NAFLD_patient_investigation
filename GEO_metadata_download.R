library(data.table)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(GEOquery)
library(rvest)
library(BiocManager)
library(limma)
library(edgeR)
library(stringi)
library(R.utils)
library(affy)
library(openxlsx)
library(WGCNA)

gds <- getGEO("GSE49541",GSEMatrix = T, getGPL = F)

metadata <- pData(gds$GSE49541_series_matrix.txt.gz) %>%
  dplyr::select(characteristics_ch1)

write.table(metadata, file = "metadata.csv")

sampleNames <- rownames(pData(gds$GSE49541_series_matrix.txt.gz))

downloader <- function(id)
{
  if(file.exists(id)){
    return(NULL)
  } else {
    getGEOSuppFiles(id)
  }
}
lapply(sampleNames, downloader)
getGEO

metadata <- read.table("~/R/tmp/NAFLD_patient_investigation/metadata.csv", quote="\"", comment.char="")


####Assemble counts####
#Use the python script "Collect_data_files" to unpack and assemble the CEl files
setwd("C:/Users/tvb217/Documents/R/tmp/NAFLD_patient_investigation/Collected_data/")
files <- list.files(pattern = ".CEL.gz", recursive = TRUE)
CEL_files <- read.affybatch(files)

eset <- rma(CEL_files)
View(eset)
geneLocat <- getBM(attributes = c("affy_hg_u133_plus_2", "external_gene_name"), 
                   values = mydata$X,
                   filters = "affy_hg_u133_plus_2",
                   mart = ensembl,
                   verbose = F)

write.exprs(eset, file = "mydata.txt")
mydata <- read.delim("~/R/tmp/NAFLD_patient_investigation/Collected_data/mydata.txt", header=T)

library(biomaRt)
deg <- AffyRNAdeg(CEL_files)
biomaRt::listFilters()
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
View(listFilters(ensembl))


View(geneLocat)

write.xlsx(geneLocat, file = "probe_to_gene_key.xlsx")



annotated_mydata <- mydata %>%
  inner_join(geneLocat, by = c("X"="affy_hg_u133_plus_2"))
View(annotated_mydata)
collapsed_mydata <- annotated_mydata
collapsed_mydata <- collapsed_mydata %>%
  distinct(X, .keep_all = T)
rownames(collapsed_mydata) <- collapsed_mydata$X
collapsed_mydata_test<-collapsed_mydata[-c(1,73)]
View(collapsed_mydata_test)
View(collapsed_mydata)


collapsed_mydata_test <- collapseRows(datET = collapsed_mydata_test, 
                    rowGroup = collapsed_mydata$external_gene_name, 
                    rowID = collapsed_mydata$X, 
                    method = "MaxMean", 
                    connectivityBasedCollapsing = F)
View(collapsed_mydata_test$datETcollapsed)



####QC of dataset####

####analysis with limma####
res <- collapsed_mydata_test$datETcollapsed
colnames(res) <- str_remove_all(colnames(res), ".CEL.gz")
View(res)
mdsData <- plotMDS(res, ndim = 3, plot = FALSE)$cmdscale.out
all(rownames(mdsData) == rownames(metadata))

#as sample GSM789122 was missing, it is removed from metadata
metadata <- metadata %>% 
  filter(rownames(metadata)!="GSM789121")
setup <- metadata %>%
  mutate(Sample = rownames(metadata))

all(colnames(res) == rownames(metadata))

mdsData <- cbind(setup,mdsData)
View(mdsData)
colnames(mdsData)[3]<-"V1"
colnames(mdsData)[4]<-"V2"
colnames(mdsData)[5]<-"V3"
colnames(mdsData)[1]<-"Group"

View(mdsData)


  ggplot(mdsData, aes(x = V1, y = V2, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")
  
  ggplot(mdsData, aes(x = V1, y = V3, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")
  
  ggplot(mdsData, aes(x = V2, y = V3, colour = Group)) +
    geom_point() +
    scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")
#limma analysis
design <- model.matrix(~ 0 + characteristics_ch1, metadata)
colnames(design)[1]<-"Stage_3_4"
colnames(design)[2]<-"Stage_0_1"


fit <- lmFit(res, design = design, method = "robust")
cont.matrix <- makeContrasts(NASH=Stage_3_4-Stage_0_1, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)  
fit2 <- eBayes(fit2)
View(fit2)

resultTable <- topTable(fit2, adjust.method = "BH", coef = 1, number = Inf)

write.xlsx(resultTable, file = "limma_results_collapsed.xlsx")  
View(resultTable)
####Go analysis####
NASH_list_significant <- resultTable %>%
  filter(adj.P.Val<0.05)
View(NASH_list_significant)

library(clusterProfiler)
library(org.Hs.eg.db)
eg= bitr(rownames(NASH_list_significant), 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = "org.Hs.eg.db",
         drop = T)

View(eg)

bg = bitr(rownames(resultTable), 
          fromType = "SYMBOL", 
          toType = "ENTREZID", 
          OrgDb = "org.Hs.eg.db",
          drop = T)
goResults <- enrichGO(gene = eg$ENTREZID,
                      universe = bg$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP")
dotplot(goResults)

#NAD 	GO:0009435
View(goResults)
cnetplot(goResults)
goResults <- setReadable(goResults, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
cnetplot(goResults)

goTerms <- goResults@result %>%
  filter(ID == "GO:0009435"	)
View(goTerms)
NADGenes <- goTerms$geneID
NADGenes_list <- unlist(strsplit(NADGenes, "/"))

significant_NAD_terms <- resultTable %>%
  filter(rownames(resultTable) %in% NADGenes_list)
View(significant_NAD_terms)

Selected_candidates <- c("NAMPT", 
                         "NMNAT1", 
                         "NMNAT2", 
                         "NMNAT3", 
                         "NMRK1", 
                         "NMRK2", 
                         "NADSYN1", 
                         "NADSYN", 
                         "TDO2", 
                         "IDO", 
                         "NAPRT", 
                         "SIRT1", 
                         "SIRT3", 
                         "PARP1", 
                         "CD38", 
                         "NNMT", 
                         "NADK", 
                         "HADH", 
                         "KMO",
                         "AFMID")
NAD_screen <- resultTable %>%
  filter(rownames(resultTable) %in% Selected_candidates)
View(NAD_screen)

ggplot(NAD_screen, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point()

