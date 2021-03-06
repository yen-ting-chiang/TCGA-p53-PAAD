#TCGAbiolinks(TCGA RNAseq data download)-----------------------------------

#BiocManager::install("TCGAbiolinks")
#browseVignettes("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("DT")

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(project = "TCGA-PAAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query,
            method = "api",
            files.per.chunk = 10)

RNAseq_data <- GDCprepare(query,
                          save = TRUE,
                          save.filename = "PAAD_RNAseq_counts.rda")


# datatable(as.data.frame(colData(data)),
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = FALSE)
# datatable(assay(data)[1:100,],
#           options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
#           rownames = TRUE)
#rowRanges(data)


#deseq2 RowData prepare----------------------------------------------------

# load("PAAD_RNAseq_counts.rda")

RNAseq_data_matrix= as.data.frame(assay(data))
save(RNAseq_data_matrix,
     file = "PAAD_RNAseq_data_matrix.rdata")
write.csv(RNAseq_data_matrix, 
          file="PAAD_RNAseq_data_matrix.csv")

# load("PAAD_RNAseq_data_matrix.rdata")

RNAseq_data_matrix_sample <- RNAseq_data_matrix
colnames(RNAseq_data_matrix_sample)=
  gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
       "",
       colnames(RNAseq_data_matrix_sample))


write.csv(RNAseq_data_matrix_sample, 
          file="RNAseq_data_matrix_sample.csv",
          quote=FALSE)



#TCGAbiolinks(TCGA mutation data download)------

PAAD.muse.maf <- GDCquery_Maf("PAAD", 
                              pipelines = "muse",
                              save.csv=TRUE)
save(PAAD.muse.maf,
     file = "PAAD.muse.maf.rdata")


load("PAAD.muse.maf.rdata")
#PAAD.muse.maf=read.csv("TCGA.PAAD.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf.csv",  header=FALSE, stringsAsFactors = FALSE)


PAAD.muse.maf_TP53_clean=PAAD.muse.maf %>% 
  filter(Hugo_Symbol=="TP53") %>% 
  mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  select(Tumor_Case_Barcode,
         Protein_position,
         HGVSp_Short,
         Exon_Number,
         Variant_Classification,
         Consequence,
         IMPACT,
         SIFT,
         PolyPhen,
  )

write.csv(PAAD.muse.maf_TP53_clean, 
          file="PAAD.muse.maf_TP53_clean.csv")

#create list of all case IDs with mutation data-------------------------------
PAAD.muse.maf_case_list=PAAD.muse.maf%>%
  mutate(Tumor_Case_Barcode=gsub("\\S-\\S\\S\\S-\\S\\S\\S\\S-\\S\\S$", 
                                 "",
                                 Tumor_Sample_Barcode))%>%
  select(Tumor_Case_Barcode)

PAAD.muse.maf_case_list_unique <- PAAD.muse.maf_case_list
PAAD.muse.maf_case_list_unique=unique(PAAD.muse.maf_case_list_unique)
write.csv(PAAD.muse.maf_case_list_unique, file = "PAAD.muse.maf_case_list_unique.csv")



#create the RowData that matches the ColData--------------------------------------------------------

#delete duplicated RNAseq data
RNAseq_data_matrix_sample_delete_duplicated <- 
  RNAseq_data_matrix_sample[, !duplicated(colnames(RNAseq_data_matrix_sample))]

mutation_list <- PAAD.muse.maf_case_list_unique[,1, drop = TRUE]

#select samples that have mutation data
RowData <- 
  RNAseq_data_matrix_sample_delete_duplicated[, colnames(RNAseq_data_matrix_sample_delete_duplicated)%in%
                                                (mutation_list)]

# save(RowData,
#      file = "RowData.rdata")
# save(PAAD.muse.maf_TP53_clean,
#      file = "PAAD.muse.maf_TP53_clean.rdata")
# save(PAAD.muse.maf_case_list_unique,
#      file = "PAAD.muse.maf_case_list_unique.rdata")
# write.csv(RowData, file = "RowData.csv")

load("RowData.rdata")
load("PAAD.muse.maf_TP53_clean.rdata")
laod("PAAD.muse.maf_case_list_unique.rdata")


#create the TP53-WT sample Coldata -----------

samples_with_TP53_MT_list <- PAAD.muse.maf_TP53_clean[,1, drop = TRUE]
samples_with_TP53_MT_list_unique <- unique(samples_with_TP53_MT_list)


DF1 <- data.frame(mutation_list)
DF1.1 <- DF1 %>% filter(mutation_list %in% 
                          samples_with_TP53_MT_list_unique == FALSE) %>% 
  mutate(Protein_position = "WT", HGVSp_Short = "WT", IMPACT = "WT")
colnames(DF1.1)[1] <- "Tumor_Case_Barcode"

#create the TP53-double-mutation sample Coldata -----------

DF2 <- data.frame(mutation_list)
DF2.1 <- PAAD.muse.maf_TP53_clean %>% 
  filter(duplicated(Tumor_Case_Barcode) == TRUE)%>% 
  select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT)%>% 
  mutate(Protein_position = "DB", HGVSp_Short = "DB", IMPACT = "DB")

#create the TP53-single-mutation sample Coldata -----------
DF3 <-  PAAD.muse.maf_TP53_clean%>% 
  select(Tumor_Case_Barcode, Protein_position, HGVSp_Short, IMPACT) %>% 
  filter(Tumor_Case_Barcode %in% DF2.1[,1, drop = TRUE] == FALSE)


#combine ColDatas----------------------------------------------------
DF4 <- full_join(DF1.1, DF2.1, 
                 by = c("Tumor_Case_Barcode", 
                        "Protein_position", 
                        "HGVSp_Short", 
                        "IMPACT"))

DF5 <- full_join(DF4, DF3, 
                 by = c("Tumor_Case_Barcode", 
                        "Protein_position", 
                        "HGVSp_Short", 
                        "IMPACT"))

ColData <- DF5 %>% 
  filter(Tumor_Case_Barcode %in% colnames(RowData) == TRUE)

# save(ColData, file = "ColData.rdata")
# write.csv(ColData, file = "ColData.csv")
load("ColData.rdata")


#launch DESeq2 analysis----------------------------------------------

load("RowData.rdata")
load("ColData.rdata")

library("DESeq2")
#cts <- RowData_match

dds <- DESeqDataSetFromMatrix(countData = as.matrix(RowData),
                              colData = ColData,
                              design= ~ IMPACT)

dds$IMPACT <- relevel(dds$IMPACT, ref = "WT")
dds <- DESeq(dds)


resultsNames(dds) # lists the coefficients
res1 <- results(dds, name="IMPACT_HIGH_vs_WT")

save(res1, file = "IMPACT_HIGH_vs_WT.rdata")
write.csv(as.data.frame(res1), 
          file="IMPACT_HIGH_vs_WT.csv")

res2 <- results(dds, name="IMPACT_MODERATE_vs_WT")

save(res2, file = "IMPACT_MODERATE_vs_WT.rdata")
write.csv(as.data.frame(res2), 
          file="IMPACT_MODERATE_vs_WT.csv")

res3 <- results(dds, name="IMPACT_MODERATE_vs_WT")

write.csv(as.data.frame(res3),
          file="IMPACT_MODERATE_vs_WT.csv")

#launch GSEA analysis----------------------------------------------------

library(dplyr)
library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)

D <- as.data.frame(res2)
head(D)
str(D)
summary(D)

names(D)[1] <- "ENSEMBL_ID"
D2=D%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)


D2_ENTREZID<- bitr(rownames(D2), fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)

D2 <- tibble::rownames_to_column(D2, "ENSEMBL")

D3=right_join(D2,D2_ENTREZID, by=c("ENSEMBL"="ENSEMBL"))
head(D3)
str(D3)
summary(D3)

## feature 1: numeric vector
geneList_ENTREZ <- D3[,3]
## feature 2: named vector
names(geneList_ENTREZ) <- D3[,8]
## feature 3: decreasing order
geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)
head(geneList_ENTREZ)
str(geneList_ENTREZ)

keytypes(org.Hs.eg.db)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

H=GSEA(
  geneList_ENTREZ,
  exponent = 1,
  minGSSize = 15,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  TERM2GENE= m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)
write.table(as.data.frame(H@result), 
            file="PAAD_IMPACT_MODERATE_vs_WT_ENTERZ_padjNAcut_H_25.txt",
            sep="\t", 
            row.names = FALSE,
            col.names = TRUE,
            quote=FALSE)


# launch DEG analysis-----------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

d <- read.csv("IMPACT_HIGH_vs_WT.csv")
head(d)
str(d)
summary(d)

names(d)[1] <- "ENSEMBL_ID"

#positive
d2=d%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange>=0)%>%
  dplyr::filter(padj<0.05)
d2_ENTREZID<- bitr(d2[,1], fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)
d2_SYMBOL<- bitr(d2[,1], fromType = "ENSEMBL",
                 toType = c("SYMBOL"),
                 OrgDb = org.Hs.eg.db,
                 drop = TRUE)
d3=left_join(d2,d2_ENTREZID, by=c("ENSEMBL_ID"="ENSEMBL"))
d4=left_join(d3,d2_SYMBOL, by=c("ENSEMBL_ID"="ENSEMBL"))
d4<- d4[ , c(1,8, 9, 2:7)]
d5=d4%>%
  arrange(desc(log2FoldChange))
write.csv(d5,file="PAAD_IMPACT_HIGH_vs_WT_DEG_pos_unfilter.csv",
          row.names = FALSE,
          quote=FALSE)

#negative
d12=d%>%
  dplyr::filter(ENSEMBL_ID !="__no_feature"&
                  ENSEMBL_ID !="__ambiguous"& 
                  ENSEMBL_ID !="__too_low_aQual"&
                  ENSEMBL_ID !="__not_aligned"&
                  ENSEMBL_ID !="__alignment_not_unique")%>%
  dplyr::filter(is.na(log2FoldChange)==FALSE)%>%
  dplyr::filter(is.na(padj)==FALSE)%>%
  dplyr::filter(log2FoldChange<=0)%>%
  dplyr::filter(padj<0.05)
d12_ENTREZID<- bitr(d12[,1], fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db,
                    drop = TRUE)
d12_SYMBOL<- bitr(d12[,1], fromType = "ENSEMBL",
                  toType = c("SYMBOL"),
                  OrgDb = org.Hs.eg.db,
                  drop = TRUE)
d13=left_join(d12,d12_ENTREZID, by=c("ENSEMBL_ID"="ENSEMBL"))
d14=left_join(d13,d12_SYMBOL, by=c("ENSEMBL_ID"="ENSEMBL"))
d14<- d14[ , c(1,8, 9, 2:7)]
d15=d14%>%
  arrange(log2FoldChange)
write.csv(d15,file="PAAD_IMPACT_HIGH_vs_WT_DEG_neg_unfilter.csv",
          row.names = FALSE,
          quote=FALSE)
