##### UBE2H in TCGA, DepMap, Yeast, Stingele, Durrbaum  ####
# UBE2H in TCGA DepMap yeast_KMS.R
# for Figure 6 and supplementary Figure 7
# analyze UBE2h in other aneuploid conditions: TCGA, DepMap, Yeast, other paired screens. etc. 

### 2025.12.12
### Author: Klaske M. Schukken

### Overview: Get TCGA Data: UBE2H expression with tissue type, 7q copy number
##            also get DepMap data. See UBE2H expression per 7q copy number. correlation with aneuploidy score
##            For Supplementary Figure 7 and parts of Figure 6


### Expression and protein data in aneuploid yeast 
## Yeast Data from: Dephoure et al. 2014, eLife; 3:e03023

## Paired human data: Stingele et al. 2012

## Durrbaum et al. 2014 Protein data 


## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"



##### Libraries #####

library(ggplot2) 
library(reshape2)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')
library(xlsx)
library("readr")
library("reshape2")
library(survival)
library(ggsurvfit)
library(dplyr)

##### Get Data #####
setwd(Dependency)

## TCGA
# 1-s2.0-S1535610818301119-mmc2.xlsx
# TCGA human tumor data, from Alison Taylor et al. (2018) 
# Genomic and Functional Approaches to Understanding Cancer Aneuploidy. 
# Volume 33, Issue 4, 9 April 2018, Pages 676-689.e3
# https://www.sciencedirect.com/science/article/pii/S1535610818301119#app2, Table S2
# I removed the top row of Table S2, which was just the name, so the actual column names would be the headers. 
# Downloaded November 10, 2023
TCGA_Aneuploidy<- read.csv("1-s2.0-S1535610818301119-mmc2.csv", header=TRUE)

# TCGA human tumor data, patient data
# Compiled by Ryan Hagenson, November 2023 
# Merges TCGA data, downloaded November 14, 2023
TCGA_Patient <- read.delim("TCGA-patients.tsv", header=TRUE, sep="\t")


# TCGA human tumor data, Sample data
# Compiled by Ryan Hagenson, November 2023
# Merges TCGA data, downloaded November 14, 2023
TCGA_Sample <- read.delim("TCGA-samples.tsv", header=TRUE, sep="\t")



# TCGA RNA expression data for UBE2H 
# Downloaded December 12, 2025  from cBioportal
# mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)
# https://www.cbioportal.org/results/download?case_set_id=all&gene_list=UBE2H&cancer_study_list=5c8a7d55e4b046111fee2296&plots_horz_selection=%7B%22dataType%22%3A%22clinical_attribute%22%2C%22selectedDataSourceOption%22%3A%22DFS_MONTHS%22%7D&plots_vert_selection=%7B%22selectedGeneOption%22%3A7328%7D&plots_coloring_selection=%7B%22colorByMutationType%22%3A%22false%22%2C%22colorByCopyNumber%22%3A%22true%22%2C%22colorBySv%22%3A%22false%22%7D
TCGA_RNA.UBE2H<- read.delim("TCGA UBE2H mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt", header= TRUE, sep="\t")


# transcription cluster 76- metabolism #
# downloaded from cBioportal HPA on January 14, 2025
# genes in cluster 76- non tissue specific metabolism 
Cluster76<- read.delim("expressionclustertissue_76_Non-specific.tsv", header= TRUE, sep="\t")




### DEPMAP ##
# Model.csv
# Sample info about DepMap cell lines
# Metadata for all of DepMap’s cancer models/cell lines. A full description of each column is available in the DepMap Release README file.
#Current DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions
#DepMap, Broad (2024). DepMap 24Q2 Public. Figshare+. Dataset. https://doi.org/10.25452/figshare.plus.25880521.v1

# Downloaded November 8, 2024
# 24Q2
Cell_Info<- read.csv("Model.csv", header=TRUE)


# aneuploidy_scores.csv
# Aneuploidy data per cell line.  Including total aneuploidy score, and cell line mean ploidy and +/- Genome doubling
# Downloaded November 9, 2023, from DepMap
# May 2023 data
Aneuploidy_Scores<- read.csv("aneuploidy_scores.csv", header=TRUE)


# arm_call_scores.csv
# depmap aneuploidy chromosome copy number dataset 
arm_call_scores<- read.csv("arm_call_scores.csv", header= TRUE)


# CRISPRGeneEffect.csv
# CRISPR screen data from DepMap
# Downloaded November 9, 2023
# May 2023 data
CRISPR_Gene_all<- read.csv("CRISPRGeneEffect.csv", header=TRUE)
colnames(CRISPR_Gene_all)[which(names(CRISPR_Gene_all) == "X")] <- "ModelID"


# OmicsExpressionProteinCodingGenesTPMLogp1.csv
# DepMap RNA expression files
# DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. Dataset.
# downloaded April 9, 2025
# DepMap Public 24Q4
DepMap_Expression<- read.csv("OmicsExpressionProteinCodingGenesTPMLogp1.csv", header=TRUE)




# DepMap Proteomics data
# downloaded August 2020
# Normalized Protein Expression data 
# https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
# Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia

Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns

m <- Protein_Expression.1$Protein_Id   
ID_Protein_Expression2<- as.data.frame(t(as.matrix(Protein_Expression.1[,-1] ))) #switch rows and collumns
colnames(ID_Protein_Expression2) <- m

ID_Protein_Expression2$Cell_Lines<-factor(rownames(ID_Protein_Expression2)) 

for (i in 1:(length(ID_Protein_Expression2)-1)){
  ID_Protein_Expression2[,i]<- as.numeric(ID_Protein_Expression2[,i])
}

ID_Protein_Expression2$Cell_Lines<-gsub("_TenPx.*","", ID_Protein_Expression2$Cell_Lines)



CCLE.Depmap.names<-read.delim2("DepMap-2018q3-celllines.csv", 
                               dec=".", header = TRUE, sep=",")
CCLE.Depmap.name<-CCLE.Depmap.names[,1:2] #get CCLE names and depmap ID correlated table

ID_Protein_effect_4<-merge(x= CCLE.Depmap.name, y= ID_Protein_Expression2, 
                           by.x="CCLE_Name", by.y="Cell_Lines", 
                           sort = TRUE)



# Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
#This data is from protein paper, not HGNC
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")




###
## Compare chromosome loss/gain difference with protein expression changes. 
## Data from Dephoure et al. 2014, eLife; 3:e03023
## Quantitative proteomic analysis reveals posttranslational responses to aneuploidy in yeast
## This is the data Jason M. Sheltzer prepared. 

# Yeast TMT means: 
Yeast.Prot<- read_xlsx("yeast-human aneuploidy.xlsx") 


##  Get Stingele et al. 2012  HCT116 aneuploidy data ##
Stingele.HCT116<- read_xls("Stingele.etal.2012.Genomics.xls", sheet=2)
Stingele.HCT116$Gene<- sub(";.*", "", Stingele.HCT116$`Gene Names`) 

Stingele<- Stingele.HCT116[,c(3,4,9,20, 26, 30, 34,45)]


### 
# Durrbaum et al. 2014 Unique features of the transcriptional response to model aneuploidy in human cells 
# BMC Genomics. 2014 Feb 18;15:139. doi: 10.1186/1471-2164-15-139
# HCT116 aneuploid clones RNA data
H.R.RNA<- read.xlsx("1471-2164-15-139-S2.xlsx", sheetIndex = 1)




#### Durrbaum et al. 2014 HCT116 RNA #####

H.R.RNA
# "hCG_18056;tcag7.352;Ube2h;UBE2H"
H.R.RNA$log2HCT116.5_4_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.5_4_2N))
H.R.RNA$p.value_HCT116.5_4<- as.numeric(as.character(H.R.RNA$p.value_HCT116.5_4))

H.R.RNA$log2HCT116.H2B.GFP.5_4_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.H2B.GFP.5_4_2N))
H.R.RNA$p.value_HCT116.H2B.5_4<- as.numeric(as.character(H.R.RNA$p.value_HCT116.H2B.5_4))

H.R.RNA$log2HCT116.3_3_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.3_3_2N))
H.R.RNA$p.value_HCT116.3_3<- as.numeric(as.character(H.R.RNA$p.value_HCT116.3_3))

H.R.RNA$log2HPT1_2N<- as.numeric(as.character(H.R.RNA$log2HPT1_2N))
H.R.RNA$p.value.HPT1<- as.numeric(as.character(H.R.RNA$p.value.HPT1))

H.R.RNA$log2HPT2_2N<- as.numeric(as.character(H.R.RNA$log2HPT2_2N))
H.R.RNA$p.value.HPT2<- as.numeric(as.character(H.R.RNA$p.value.HPT2))

H.R.RNA$log2.RPE1.H2B.GFP.21_3_2N<- as.numeric(as.character(H.R.RNA$log2.RPE1.H2B.GFP.21_3_2N))
H.R.RNA$p.value_RPE1.H2B.21_3<- as.numeric(as.character(H.R.RNA$p.value_RPE1.H2B.21_3))

H.R.RNA$log2.RPE1.5_3.12_3_2N<- as.numeric(as.character(H.R.RNA$log2.RPE1.5_3.12_3_2N))
H.R.RNA$p.value_RPE1.12_3.5_3<- as.numeric(as.character(H.R.RNA$p.value_RPE1.12_3.5_3))




H.R.RNA$meanlog2FC_H<- rowMeans(H.R.RNA[, c(3,6,9)] )
H.R.RNA$RankMeanFC_H<- rank(H.R.RNA$meanlog2FC_H)

ggplot(H.R.RNA, aes(x= RankMeanFC_H, y= meanlog2FC_H))+
  geom_point()+
  xlim(0,13847)+
  xlab("Rank HCT116 Aneu (Durrbaum et al.)")+
  ylab("log2FC RNA")+
  geom_point(data= subset(H.R.RNA, Gene.name == "hCG_18056;tcag7.352;Ube2h;UBE2H"), color="red", size=3)+
  theme_classic()
# Durrbaum et al. 2014 
# plot.Durrbaum.RNA.AneuFC.rank.UBE2H_HCT116Ts5.Ts3.Tetra5.pdf
# Figure S7A
# UBE2H is in top 7.2% of aneuploid dependencies in HCT116. 



#### Stingele et al. 2012 HCT116 protein ####

colnames(Stingele)

Stingele[,3]<- as.numeric(Stingele$`mRNA HCT116 5/4`)
Stingele[,4]<- as.numeric(Stingele$`Protein HCT116 5/4`)
Stingele[,5]<- as.numeric(Stingele$`Protein HCT116 H2B-GFP 5/4`)
Stingele[,6]<- as.numeric(Stingele$`Protein HCT116 H2B-GFP 5/3`)
Stingele[,7]<- as.numeric(Stingele$`Protein HCT116 3/3`)

Stingele$MeanHCT116<- rowMeans(Stingele[,c(4,5,6,7)])

Stingele$RankH<- rank(Stingele$MeanHCT116)

ggplot(Stingele, aes(x= RankH, y= MeanHCT116))+
  geom_point()+
  xlim(0,3564)+
  xlab("Rank HCT116 Aneu (Stingele et al.)")+
  ylab("log2FC Protein")+
  geom_point(data= subset(Stingele, Gene == "UBE2H"), color="red", size=3)+
  theme_classic()
# Plot.Stingele.HCT116.Protein.Rank.pdf 
# UBE2H is 3552 out of 3564 top 99.6% of upregulated proteins in aneuploid HCT116

# Figure S7B



#### Yeast: Protein & RNA difference in aneuploid yeast  #####

###Plot difference upon trisomy in yeast 
subset(Yeast.Prot, `Yeast Gene`== "UBC8")

ggplot(Yeast.Prot, aes(x=`Yeast RNA`, y=`Yeast TMT - protein`))+
  geom_point(color="black")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red", size=2)+
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  ylab("Protein difference (Dephoure et al.)\nin aneuploid yeast")+
  xlab("RNA difference (Dephoure et al.)\nin aneuploid yeast")
# size 4x4 
# plot.Dephoure.YeastTrisomyExpression_UBE2H


Yeast.Prot$RankYeastRNA<- rank(Yeast.Prot$`Yeast RNA`)
ggplot(Yeast.Prot, aes(x= RankYeastRNA, y= `Yeast RNA`))+
  geom_point()+
  xlab("Rank Yeast Aneuploid (Dephoure et al.)")+
  ylab("log2FC RNA")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red")+
  theme_classic()
# size 4x4 
# Figure S7D
# plot.Dephoure.YeastTrisomyRNA_UBE2H_Rank
# 665/738= top 9.9% of RNA


Yeast.Prot$RankYeast<- rank(Yeast.Prot$`Yeast TMT - protein`)
ggplot(Yeast.Prot, aes(x= RankYeast, y= `Yeast TMT - protein`))+
  geom_point()+
  xlab("Rank Yeast Aneuploid (Dephoure et al.)")+
  ylab("log2FC Protein")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red")+
  theme_classic()
# size 4x4 
# plot.Dephoure.YeastTrisomyProtein_UBE2H_Rank
# 700/738= top 5.2% of proteins
# Figure S7C





#### TCGA UBE2H RNA expression #####
# Merge data to get more info and about patients
TCGA_Patient_info<- merge(TCGA_RNA.UBE2H, TCGA_Sample[, c(1,2,4, 13)], 
                                by.x= "SAMPLE_ID", by.y= "SAMPLE_ID")
TCGA_Patient_info<- merge(x= TCGA_Patient_info, y=TCGA_Patient, 
                                by.y= "PATIENT_ID", by.x= "PATIENT_ID")


ggplot(TCGA_Patient_info, aes(y=log2(UBE2H), x=CANCER_TYPE_ACRONYM))+
  geom_boxplot(outlier.shape=NA)+  
  #ylim(0,5000)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_classic()
# TCGA.UBE2H.RNA.perTissue.pdf

#Note Leukemias and Mature B-cell Neoplasms have much lower UBE2H expression than other cell lines 
# DLBC and LAML


# Aneuploidy correlation:
TCGA_Patient_info_Aneu<- merge(x= TCGA_Patient_info, y=TCGA_Aneuploidy, 
                          by.x= "SAMPLE_ID", by.y= "Sample")
TCGA_Patient_info_Aneu$PloidEst<- TCGA_Patient_info_Aneu$Genome_doublings
TCGA_Patient_info_Aneu <- TCGA_Patient_info_Aneu %>%
  mutate(PloidEst = recode(PloidEst, # make estimated ploidy: 0 doublings = 2, 1 doubling = tetraploid 4, 2 doubling is octoploid
                    `0` = 2,
                    `1` = 4,
                    `2` = 8))
TCGA_Patient_info_Aneu$AneuScore.ploidy<- TCGA_Patient_info_Aneu$AneuploidyScore.AS. / TCGA_Patient_info_Aneu$PloidEst

ggplot(TCGA_Patient_info_Aneu, aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy.ploidy.pdf
# Figure 6L

x<- cor.test(log2(TCGA_Patient_info_Aneu$UBE2H), TCGA_Patient_info_Aneu$AneuScore.ploidy)
x$p.value
x
# cor= 0.1028087
# p= 4.149232e-24


# UBE2H expression by 7q status 
ggplot(data = subset(TCGA_Patient_info_Aneu, ! X7q == "NA"), # remove samples where 7q CN status is NA
       aes(y=log2(UBE2H), x=as.factor(X7q)))+
  geom_boxplot()+
  xlab("7q Copy Number")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.7qCNboxplot.pdf
# Figure S7H
mean(subset(TCGA_Patient_info_Aneu, X7q == "0" & ! UBE2H %in% c(NA))$UBE2H) - mean(subset(TCGA_Patient_info_Aneu, X7q == "1" & ! UBE2H %in% c(NA))$UBE2H)
x<- t.test(subset(TCGA_Patient_info_Aneu, X7q == "0")$UBE2H, subset(TCGA_Patient_info_Aneu, X7q == "1")$UBE2H)
x$p.value # neutral gain p= 1.894639e-39
mean(subset(TCGA_Patient_info_Aneu, X7q == "0" & ! UBE2H %in% c(NA))$UBE2H) - mean(subset(TCGA_Patient_info_Aneu, X7q == "-1" & ! UBE2H %in% c(NA))$UBE2H)
x<- t.test(subset(TCGA_Patient_info_Aneu, X7q == "0")$UBE2H, subset(TCGA_Patient_info_Aneu, X7q == "-1")$UBE2H)
x$p.value # neutral loss 5.852346e-09


# correlation of UBE2H expression and aneuploidy score in 7q neutral cells only
ggplot(subset(TCGA_Patient_info_Aneu, X7q==0), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score (neutral 7q only)")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy.ploidy_Neutral7q.pdf
# Figure S7J
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, X7q==0)$UBE2H), subset(TCGA_Patient_info_Aneu, X7q==0)$AneuScore.ploidy)
x$p.value
x
# p=5.311165e-20 
# corr= 0.1188832




# TCGA  UBE2H RNA correspond with aneuscore per tissue
UBE2H_RNA.Aneu_Tissue<- data.frame (TissueType = character(), 
                                     UBE2H.AnCorr = numeric(), 
                                    UBE2H.An.p = character())

for (i in 1: length(unique(TCGA_Patient_info_Aneu$CANCER_TYPE_ACRONYM))) {
  Tissue<- unique(TCGA_Patient_info_Aneu$CANCER_TYPE_ACRONYM)[i]
  Corr<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == Tissue)$UBE2H), 
                  subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == Tissue)$AneuScore.ploidy, 
                  method = "pearson")
    
  UBE2H_RNA.Aneu_Tissue<- rbind(UBE2H_RNA.Aneu_Tissue, data.frame(TissueType = Tissue, 
                                                                  UBE2H.AnCorr = Corr$estimate, 
                                                                  UBE2H.An.p = Corr$p.value ))
  
}
UBE2H_RNA.Aneu_Tissue
# Most significant tissue types: 

ggplot(UBE2H_RNA.Aneu_Tissue, aes(x=UBE2H.AnCorr, y= -log2(UBE2H.An.p), label=TissueType))+
  geom_point()+
  geom_hline(yintercept=-log2(0.05))+
  geom_text(hjust=0, vjust=0)+
  theme_classic()
# plot.TCGA.TissueType.corr.pvalue.pdf
# Figure S7K


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_KIRC.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC")$AneuScore.ploidy)
x
# 0.2447174
# p= 5.45e-08
# Figure S7L

ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_ACC.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC")$AneuScore.ploidy, method = "pearson")
x
# 0.4418611
# p= 6.444e-05
# Figure 6M


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_Pancreatic.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD")$AneuScore.ploidy)
x
# 0.2145797
# p= 0.006781


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_SKCM.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM")$AneuScore.ploidy)
x
# 0.187036
# p= 9.028e-05



ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM %in% c("READ")), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_Pancreatic.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "READ")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "READ")$AneuScore.ploidy)
x
# 0.2272697
# p= 0.004589






  ## TCGA Genes in cluster 76 Metabolism  #####


MitoTranslationTranscription <- c("TRMT10C","HSD17B10","PRORP", 
                                  "POLRMT", "TFAM", "TFB2M",
                                  "MRPL58", "MTRFR", "MTRF1", "MTRF1L", 
                                  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                  "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", 
                                  "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", 
                                  "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", 
                                  "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", 
                                  "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B", 
                                  "FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "TBRG4", "FASTKD5")
MitoElectronChainnAssembly <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2", "NDUFAB1","NDUFA1",
                                "NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
                                "NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                                "NDUFB11","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3", 
                                "SDHA", "SDHB", "SDHC", "SDHD", 
                                "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB","CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11", 
                                "COX4I1","COX4I2","COX5A","COX5B","COX6A1", "COX6A2","COX6B1","COX6B2","COX6C", "COX7A1","COX7A2","COX7B","COX7B2","COX7C",
                                "COX8A","COX8C","MT-CO1","MT-CO2","MT-CO3", 
                                "ATP5F1A", "ATP5F1B", "ATP5F1C",  "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MJ", "ATP5MK", 
                                "MT-ATP6", "MT-ATP8", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5IF1", 
                                "ATPAF1","ATPAF2","BCS1L","CHCHD7","CMC1","CMC2", "COA1",  "COA3", "COA4","COA5", "COA6",  "COA7", "COA8",
                                "COX7A2L","COX10", "COX11", "COX14","COX15", "COX16",  "COX17", "COX18", "COX19", "COX20",
                                "DMAC1", "DMAC2", "FDXR", "FDX2", "FMC1", "FOXRED1","HCCS","HIGD1A","HIGD2A","LRPPRC", "LYRM7",
                                "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6","NDUFAF7",  "NDUFAF8","NUBPL",
                                "OCIAD2", "OXA1L", "PET100","PET117","PNKD", "SCO1", "SCO2",  "SDHAF1","SDHAF2","SDHAF3", "SDHAF4", "SFXN4",
                                "SMIM20", "SURF1","TACO1", "TIMMDC1","TMEM70","TMEM126A", "TMEM177","TMEM186","TMEM223","TMEM242","TTC19", 
                                "UQCC1","UQCC2","UQCC3", "UQCC4","UQCC5","UQCC6", 
                                "ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186", 
                                "BOLA3", "FDXR", "FXN", "GLRX5", "IBA57", "ISCA1", "ISCA2", "ISCU", "LYRM4", "NFS1", "NFU1", "NUBPL") # HGNC groups

length(Cluster76$Gene[Cluster76$Gene %in% MitoElectronChainnAssembly]) #14 genes
length(Cluster76$Gene[Cluster76$Gene %in% MitoTranslationTranscription]) #17


  ## TCGA Percent 7q gain #####
# Merge data to get more info and about patients
TCGA_Patient_info<- merge(TCGA_RNA.UBE2H, TCGA_Sample[, c(1,2,4, 13)], 
                          by.x= "SAMPLE_ID", by.y= "SAMPLE_ID")
TCGA_Patient_info<- merge(x= TCGA_Patient_info, y=TCGA_Patient, 
                          by.y= "PATIENT_ID", by.x= "PATIENT_ID")
TCGA_Patient_info_Aneu<- merge(x= TCGA_Patient_info, y=TCGA_Aneuploidy, 
                               by.x= "SAMPLE_ID", by.y= "Sample")


table(TCGA_Aneuploidy$X7q)
#-1    0    1 
#507 6171 2490

# Percent of TCGA tumors with 7q gain: 
PercentGainLoss<- data.frame(Arm= character(),
                             percentGain= as.numeric(), 
                             percentLoss = as.numeric())
nTumor<- length(TCGA_Aneuploidy$X1p)

for (i in 14:52){
  TableCN <- table(TCGA_Aneuploidy[i])
  PercentGainLoss<- rbind(PercentGainLoss, data.frame(Arm= colnames(TCGA_Aneuploidy)[i],
                                                      percentGain = TableCN["1"]/nTumor, 
                                                      percentLoss= TableCN["-1"]/nTumor))
}
PercentGainLoss<- PercentGainLoss[order(PercentGainLoss$percentGain),]
PercentGainLoss
# 7q percent gain: 0.23664703  
# 7q percent loss: 0.04818476

# 7q gain 5th most common in TCGA

#### DepMap data: plot UBE2H RNA expression #####

Cell_Info

Aneuploidy_Scores[1:10,]

CRISPR_Gene_all

DepMap_Expression

Aneuploidy_Scores$Aneu.ploidy<- Aneuploidy_Scores$Aneuploidy.score/Aneuploidy_Scores$Ploidy
DepMap_RNA_Aneu<- merge(Aneuploidy_Scores, DepMap_Expression, by.x="DepMap_ID", by.y="X")
DepMap_RNA_Aneu<- merge(DepMap_RNA_Aneu, arm_call_scores, by.x="DepMap_ID", by.y= "X")
DepMap_RNA_Aneu<- merge(DepMap_RNA_Aneu, CRISPR_Gene_all, by.x="DepMap_ID", by.y="ModelID")
# note: gene ID with ".x" is expression, gene id with ".y" is CRISPR effect score

#UBE2H correlation with aneuploidy score
ggplot(DepMap_RNA_Aneu, aes(y=log2(UBE2H..7328..x), x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("DepMap: Aneuploidy Score")+ 
  ylab("UBE2H expression")+
  theme_classic()
# DepMap.UBE2H.RNA.Aneuploidy.ploidy.pdf
# Figure S7E
cor.test(DepMap_RNA_Aneu$UBE2H..7328..x, DepMap_RNA_Aneu$Aneu.ploidy)
# 0.000671
# 0.1348323



# UBE2H expression by 7q status 
ggplot(data = DepMap_RNA_Aneu, # remove samples where 7q CN status is NA
       aes(y=log2(UBE2H..7328..x), x=as.factor(X7q)))+
  geom_boxplot()+
  xlab("7q Copy Number")+ 
  ylab("UBE2H expression")+
  theme_classic()
# DepMap.UBE2H.RNA.7qCNboxplot.pdf
# Figure S7F
mean(subset(DepMap_RNA_Aneu, X7q == "0" )$UBE2H..7328..x) - mean(subset(DepMap_RNA_Aneu, X7q == "1" )$UBE2H..7328..x)
x<- t.test(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x, subset(DepMap_RNA_Aneu, X7q == "1")$UBE2H..7328..x)
x$p.value # neutral gain p= 1.718648e-12
mean(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x) - mean(subset(DepMap_RNA_Aneu, X7q == "-1" )$UBE2H..7328..x)
x<- t.test(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x, subset(DepMap_RNA_Aneu, X7q == "-1")$UBE2H..7328..x)
x$p.value # neutral loss 0.03292398
 

# correlation of UBE2H expression and aneuploidy score in 7q neutral cells only
ggplot(subset(DepMap_RNA_Aneu, X7q==0), aes(y=log2(UBE2H..7328..x), x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  ylab("UBE2H Expression")+ 
  xlab("DepMap: Aneuploidy Score")+
  theme_classic()
# DepMap.UBE2H.RNA.Aneuploidy.ploidy_Neutral7q.pdf
# Figure S7I
cor.test(log2(subset(DepMap_RNA_Aneu, X7q==0)$UBE2H..7328..x), subset(DepMap_RNA_Aneu, X7q==0)$Aneu.ploidy)
# p=0.01658 
# corr= 0.1322125 




#### DepMap data: plot UBE2H Protein #####

ID_Protein_effect_4
Cell_Info[c(1,6),]
Aneuploidy_Scores[1:10,]


Aneuploidy_Scores$Aneu.ploidy<- Aneuploidy_Scores$Aneuploidy.score/Aneuploidy_Scores$Ploidy
DepMap_Protein_Aneu<- merge(Aneuploidy_Scores, ID_Protein_effect_4, by.x="DepMap_ID", by.y="Broad_ID")
DepMap_Protein_Aneu<- merge(DepMap_Protein_Aneu, Cell_Info[, c(1,6)], by.x="DepMap_ID", by.y="ModelID")

DepMap_Protein_Aneu2<- merge(DepMap_Protein_Aneu, arm_call_scores, by.x="DepMap_ID", by.y= "X")
DepMap_Protein_Aneu2<- merge(DepMap_Protein_Aneu2, CRISPR_Gene_all, by.x="DepMap_ID", by.y="ModelID")

subset(Protein_ProID, Gene_Symbol == "UBE2H")

#UBE2H correlation with aneuploidy score
ggplot(DepMap_Protein_Aneu, aes(y=((`sp|P62256|UBE2H_HUMAN`)), x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("DepMap: Aneuploidy Score")+ 
  ylab("UBE2H protein")+
  theme_classic()
# DepMap.UBE2H.protein.Aneuploidy.ploidy.pdf
cor.test(DepMap_Protein_Aneu$`sp|P62256|UBE2H_HUMAN`, DepMap_Protein_Aneu$Aneu.ploidy)
# p-value= 0.4702
# cor= 0.03745756


# UBE2H protein by lineage
ggplot(data = DepMap_Protein_Aneu, # remove samples where 7q CN status is NA
       aes(y=(`sp|P62256|UBE2H_HUMAN`), x=as.factor(OncotreeLineage)))+
  geom_boxplot()+
  xlab("Lineage")+ 
  ylab("UBE2H protein")+
  theme_classic()
# DepMap.UBE2H.RNA.7qCNboxplot.pdf


# UBE2H protein by 7q status 
ggplot(data = DepMap_Protein_Aneu2, # remove samples where 7q CN status is NA
       aes(y=(`sp|P62256|UBE2H_HUMAN`), x=as.factor(X7q)))+
  geom_boxplot()+
  xlab("7q Copy Number")+ 
  ylab("UBE2H protein")+
  theme_classic()
# DepMap.UBE2H.Protein.7qCNboxplot.pdf
# Figure S7G

mean(subset(DepMap_Protein_Aneu2, X7q == "0" )$`sp|P62256|UBE2H_HUMAN`)  - mean(subset(DepMap_Protein_Aneu2, X7q == "1" )$`sp|P62256|UBE2H_HUMAN`) 
x<- t.test(subset(DepMap_Protein_Aneu2, X7q == "0")$`sp|P62256|UBE2H_HUMAN`, subset(DepMap_Protein_Aneu2, X7q == "1")$`sp|P62256|UBE2H_HUMAN`)
x$p.value # neutral gain p= 8.005047e-07    diff= 0.4659726
mean(subset(DepMap_Protein_Aneu2, X7q == "0")$`sp|P62256|UBE2H_HUMAN`) - mean(subset(DepMap_Protein_Aneu2, X7q == "-1" )$`sp|P62256|UBE2H_HUMAN`)
x<- t.test(subset(DepMap_Protein_Aneu2, X7q == "0")$`sp|P62256|UBE2H_HUMAN`, subset(DepMap_Protein_Aneu2, X7q == "-1")$`sp|P62256|UBE2H_HUMAN`)
x$p.value # neutral loss 0.05107126 NS  diff= 0.2640408


# Correlation of UBE2H expression and aneuploidy score in 7q neutral cells only
ggplot(subset(DepMap_Protein_Aneu2, X7q==0), aes(y=`sp|P62256|UBE2H_HUMAN`, x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  ylab("UBE2H Protein")+ 
  xlab("DepMap: Aneuploidy Score")+
  theme_classic()
# DepMap.UBE2H.RNA.Aneuploidy.ploidy_Neutral7q.pdf
# Figure S7I
cor.test(subset(DepMap_Protein_Aneu2, X7q==0)$`sp|P62256|UBE2H_HUMAN`, subset(DepMap_Protein_Aneu2, X7q==0)$Aneu.ploidy)
# p=0.2358  NS
# corr= 0.09280467 



UBE2H.Corr.Prot<- data.frame(Protein_Id= character(), 
                                   Corr= numeric(), 
                                   pvalue= numeric())

for (i in 8:12762){ 
  x<- cor.test(DepMap_Protein_Aneu$`sp|P62256|UBE2H_HUMAN`, DepMap_Protein_Aneu[,i])
  UBE2H.Corr.Prot<- rbind(UBE2H.Corr.Prot, data.frame(Protein_Id= colnames(DepMap_Protein_Aneu)[i], 
                                                                  Corr= x$estimate, 
                                                                  pvalue= x$p.value ))
}

subset(UBE2H.Corr.Prot, Protein_Id == "sp|P62256|UBE2H_HUMAN")

UBE2H.Corr.Prot<- UBE2H.Corr.Prot[-11665,] # Remove UBE2H no need to correlate UBE2H with itself

UBE2H.Corr.Prot$Protein_Id2<- gsub("\\|",".", UBE2H.Corr.Prot$Protein_Id)
UBE2H.Corr.Prot$Protein_Id2[1:10]
UBE2H.Corr.Prot<- merge(Protein_ProID, UBE2H.Corr.Prot, by.x="Protein_Id", by.y="Protein_Id2")


# plot genes that correlate with UBE2H 
ggplot(subset(UBE2H.Corr.Prot, ! Gene_Symbol %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription, RNA_Processing_rRNA, Ribosomal_Genes_noMT$Approved.symbol, Proteosome_Genes$Approved.symbol, RNA_Processing_Spliceosome)), aes(x=Corr))+
  geom_density(size=2)+
  geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% MitoElectronChainnAssembly), color = "#F8766D", size=2) + 
  geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% MitoTranslationTranscription ), color="gold2", size=2) +
  #geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% RNA_Processing_rRNA ), color="deepskyblue4", size=2) + #00BFC4
  #geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF", size=2) + #C77CFF
  geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1", size=2) + #7CAE00
  #geom_density(data = subset(UBE2H.Corr.Prot, Gene_Symbol %in% RNA_Processing_Spliceosome ), color="deepskyblue1", size=2) + #00BFC4
  xlab("Correlation (DepMap):\n UBE2H Protein & Protein X")+ 
  theme_classic()
# DepMap.UBE2H.protein.Corr.Protein.pdf





UBE2H.Corr.CRISPR<- data.frame(Protein_Id= character(), 
                               Corr= numeric(), 
                               pvalue= numeric())
for (i in 12802:30732){ 
  x<- cor.test(DepMap_Protein_Aneu$`sp|P62256|UBE2H_HUMAN`, DepMap_Protein_Aneu[,i])
  UBE2H.Corr.CRISPR<- rbind(UBE2H.Corr.CRISPR, data.frame(Protein_Id= colnames(DepMap_Protein_Aneu)[i], 
                                                          Corr= x$estimate, 
                                                          pvalue= x$p.value ))
}

UBE2H.Corr.CRISPR$Protein_Id2<- gsub("[..].*","", UBE2H.Corr.CRISPR$Protein_Id)
UBE2H.Corr.CRISPR$Protein_Id2[1:10]


# Plot genes that correlate with UBE2H 
ggplot(data=subset(UBE2H.Corr.CRISPR, ! Protein_Id2 %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription)), aes(x=Corr))+
  geom_density(size=2)+
  geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% MitoElectronChainnAssembly), color = "#F8766D", size=2) + 
  geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% MitoTranslationTranscription ), color="gold2", size=2) +
  #geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% RNA_Processing_rRNA ), color="deepskyblue4", size=2) + #00BFC4
  #geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF", size=2) + #C77CFF
  geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1", size=2) + #7CAE00
  #geom_density(data = subset(UBE2H.Corr.CRISPR, Protein_Id2 %in% RNA_Processing_Spliceosome ), color="deepskyblue1", size=2) + #00BFC4
  xlab("Correlation (DepMap):\n CRISPR Effect Score & UBE2H protein abundance")+ 
  theme_classic()
# DepMap.UBE2H.protein.Corr.CRISPR.pdf



subset(UBE2H.Corr.CRISPR, Corr< -0.2 & pvalue<0.05)$Protein_Id2

